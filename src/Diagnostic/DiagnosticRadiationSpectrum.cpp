#include "PyTools.h"
#include <iomanip>

#include "DiagnosticRadiationSpectrum.h"
#include "HistogramFactory.h"
#include "RadiationTools.h"
#include "RadiationTables.h"


using namespace std;


// Constructor
DiagnosticRadiationSpectrum::DiagnosticRadiationSpectrum(
    Params &params,
    SmileiMPI *smpi,
    Patch *patch,
    RadiationTables *radiation_tables_,
    int diagId
) : DiagnosticParticleBinningBase( params, smpi, patch, diagId, "RadiationSpectrum", false, PyUnicode_FromString( "" ), excludedAxes() )
{
    
    ostringstream name( "" );
    name << "DiagRadiationSpectrum #" << diagId;
    string errorPrefix = name.str();
    
    // Check that reference_angular_frequency_SI is correctly defined
    if (params.reference_angular_frequency_SI<=0.) {
        ERROR("DiagRadiationSpectrum requires 'reference_angular_frequency_SI' to be defined.");
    }
    
    // minimum chi beyond which the radiation spectrum is computed (uses minimum_chi_continuous)
    minimum_chi_continuous_ = radiation_tables_->getMinimumChiContinuous();
    
    // Normalization parameters
    two_third = 2./3.;
    double squared_fine_structure_constant = 5.325135447834466e-5;
    double normalized_classical_electron_time = 9.399637140638142e-24*params.reference_angular_frequency_SI;
    double factor  = two_third*squared_fine_structure_constant/normalized_classical_electron_time;
    factor *= sqrt(3.)/2./M_PI;
    
    // Histogram axes that should not be allowed
    vector<string> excluded_axes( 0 );
    excluded_axes.push_back( "a" );
    excluded_axes.push_back( "b" );
    excluded_axes.push_back( "theta" );
    excluded_axes.push_back( "phi" );
    
    // Get the "photon energy axis"
    PyObject* photon_energy_axis = PyTools::extract_py( "photon_energy_axis", "DiagRadiationSpectrum", diagId );
    ostringstream t("");
    t << errorPrefix << "photon_energy_axis : ";
    photon_axis = HistogramFactory::createAxis( photon_energy_axis, params, species, patch, excluded_axes, t.str(), false );
    total_axes++;
    dims.push_back( photon_axis->nbins );
    
    // construct the list of photon_energies
    photon_energies.resize( photon_axis->nbins );
    delta_energies.resize( photon_axis->nbins );
    double emin = photon_axis->actual_min;
    double emax = photon_axis->actual_max;
    double spacing = (emax-emin) / photon_axis->nbins;
    for( int i=0; i<photon_axis->nbins; i++ ) {
        photon_energies[i] = emin + (i+0.5)*spacing;
        if( photon_axis->logscale ) {
            photon_energies[i] = pow(10., photon_energies[i]);
            delta_energies[i] = pow(10., emin+i*spacing) * ( pow(10., spacing) - 1. );
        } else {
            delta_energies[i] = spacing;
        }
        delta_energies[i] *= factor;
    }
    
    // Calculate the size of the output array
    uint64_t total_size = (uint64_t)output_size * photon_axis->nbins;
    if( total_size > 2147483648 ) { // 2^31
        ERROR( errorPrefix << ": too many points (" << total_size << " > 2^312)" );
    }
    output_size = ( unsigned int ) total_size;
    
    // Output info on diagnostics
    if( smpi->isMaster() ) {
        MESSAGE( 2, photon_axis->info( "photon energy" ) );
    }

} // END DiagnosticRadiationSpectrum::DiagnosticRadiationSpectrum


DiagnosticRadiationSpectrum::~DiagnosticRadiationSpectrum()
{
    delete photon_axis;
} // END DiagnosticRadiationSpectrum::~DiagnosticRadiationSpectrum


// Called only by patch master of process master
void DiagnosticRadiationSpectrum::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if( !smpi->isMaster() || fileId_>0 ) {
        return;
    }
    
    if( newfile ) {
        DiagnosticParticleBinningBase::openFile( params, smpi, newfile );
        
        // write photon_energy_axis
        string str1 = "photon_energy_axis";
        ostringstream mystream( "" );
        mystream << photon_axis->min << " " << photon_axis->max << " "
                 << photon_axis->nbins << " " << photon_axis->logscale << " " << photon_axis->edge_inclusive;
        string str2 = mystream.str();
        H5::attr( fileId_, str1, str2 );
        
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
}

// run one particle binning diagnostic
void DiagnosticRadiationSpectrum::run( Patch* patch, int timestep, SimWindow* simWindow )
{

    vector<int> int_buffer;
    vector<double> double_buffer;
    
//    // Update spatial_min and spatial_max if needed
//    if( simWindow ) {
//        bool did_move = false;
//        for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
//            if( histogram->axes[i]->type == "moving_x" ) {
//                spatial_min[0] = histogram->axes[i]->min + simWindow->getXmoved();
//                spatial_max[0] = histogram->axes[i]->max + simWindow->getXmoved();
//                did_move = true;
//            }
//        }
//        if( ! did_move ) {
//            spatial_max[0] += simWindow->getXmoved() - spatial_min[0];
//            spatial_min[0] = simWindow->getXmoved();
//        }
//    }
    
    // loop species
    for( unsigned int ispec=0 ; ispec < species.size() ; ispec++ ) {
        
        Species *s = patch->vecSpecies[species[ispec]];
        unsigned int npart = s->particles->size();
        int_buffer   .resize( npart );
        double_buffer.resize( npart );
        
        fill(int_buffer.begin(), int_buffer.end(), 0);
        
        histogram->digitize( s, double_buffer, int_buffer, simWindow );
        
        // Sum the data into the data_sum
        // ------------------------------
        int ind;
        
        double gamma_inv, gamma, chi, xi, zeta, nu, cst;
        double two_third_ov_chi, increment0, increment;
        int iphoton_energy_max;
        
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            ind = int_buffer[ipart];
            if( ind<0 ) continue; // skip already discarded particles
            ind *= photon_axis->nbins;
            
            // Get the quantum parameter
            chi = s->particles->chi( ipart );
            
            // Update the spectrum only if the quantum parameter is sufficiently high
            if( chi <= minimum_chi_continuous_ ) continue;
            
            // Emitting particle energy (maximum of the spectrum)
            gamma = s->particles->LorentzFactor( ipart );
            gamma_inv = 1./gamma;
            two_third_ov_chi = two_third/chi;
            increment0 = gamma_inv * s->particles->weight( ipart );
            
            // Compute the maximum iteration of the loop on bins
            // ensures that xi<1;
            // that is no radiation corresponds to photon energy larger than the radiating particle energy
            if( photon_axis->logscale ) {
                gamma = log10( gamma );
            }
            iphoton_energy_max = int( (gamma - photon_axis->actual_min) * photon_axis->coeff );
            //iphoton_energy_max can not be greater than photon_energy_nbins
            iphoton_energy_max = min( iphoton_energy_max, photon_axis->nbins );
            
            // Loop on bins
            for( int i=0; i<iphoton_energy_max; i++ ) {
                xi   = photon_energies[i] * gamma_inv;
                zeta = xi / (1.-xi); // xi<1 is ensured above
                nu   = two_third_ov_chi * zeta;
                cst  = xi * zeta;
                increment = increment0 * delta_energies[i] * xi * RadiationTools::computeBesselPartsRadiatedPower(nu,cst);
                #pragma omp atomic
                data_sum[ind+i] += increment;
            }
            
        }
    
    }
    
} // END run



