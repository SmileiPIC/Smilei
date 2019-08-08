#include "PyTools.h"
#include <iomanip>

#include "DiagnosticRadiationSpectrum.h"
#include "HistogramFactory.h"
#include "RadiationTables.h"
#include "RadiationTools.h"


using namespace std;


// Constructor
DiagnosticRadiationSpectrum::DiagnosticRadiationSpectrum( Params &params, SmileiMPI* smpi, Patch* patch, int diagId )
{
    fileId_ = 0;

    int n_diag_rad_spectrum = diagId;

    // Check that reference_angular_frequency_SI is correctly defined
    if (params.reference_angular_frequency_SI<=0.) ERROR("DiagnosticRadiationSpectrum requires 'reference_angular_frequency_SI' to be defined.");

    // Normalization parameters
    two_third = 2./3.;
    double squared_fine_structure_constant = 5.325135447834466e-5;
    double normalized_classical_electron_time = 9.399637140638142e-24*params.reference_angular_frequency_SI;
    factor  = two_third*squared_fine_structure_constant/normalized_classical_electron_time;//*MG*params.timestep;
    factor *= 0.2756644477108960; //sqrt(3.)/2./M_PI;


    ostringstream name("");
    name << "Diagnotic Radiation Spectrum #" << n_diag_rad_spectrum;
    string errorPrefix = name.str();

    // get parameter "every" which describes a timestep selection
    timeSelection = new TimeSelection(
        PyTools::extract_py("every", "DiagRadiationSpectrum", n_diag_rad_spectrum),
        name.str()
    );

    // get parameter "flush_every" which describes a timestep selection for flushing the file
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py("flush_every", "DiagRadiationSpectrum", n_diag_rad_spectrum),
        name.str()
    );

    // get parameter "time_average" that determines the number of timestep to average the outputs
    time_average = 1;
    PyTools::extract("time_average",time_average,"DiagRadiationSpectrum",n_diag_rad_spectrum);
    if ( time_average < 1 ) time_average=1;
    if ( time_average > timeSelection->smallestInterval() )
        ERROR(errorPrefix << ": `time_average` is incompatible with `every`");

    // get parameter "species" that determines the species to use (can be a list of species)
    vector<string> species_names;
    if (!PyTools::extract("species",species_names,"DiagRadiationSpectrum",n_diag_rad_spectrum))
        ERROR(errorPrefix << ": parameter `species` required");
    // verify that the species exist, remove duplicates and sort by number
    species = params.FindSpecies(patch->vecSpecies, species_names);


    // start reading the "photon energy axis"
    // --------------------------------------
    PyObject* photon_energy_axis = PyTools::extract_py("photon_energy_axis", "DiagRadiationSpectrum", n_diag_rad_spectrum);
    // Axis must be a list
    if (!PyTuple_Check(photon_energy_axis) && !PyList_Check(photon_energy_axis))
        ERROR(errorPrefix << ": photon_energy_axis must be a list");
    PyObject* seq = PySequence_Fast(photon_energy_axis, "expected a sequence");

    // Axis must have 3 elements or more
    unsigned int lenAxisArgs=PySequence_Size(seq);
    if (lenAxisArgs<3)
        ERROR(errorPrefix << ": photon_energy_axis must contain at least 3 arguments (contains only " << lenAxisArgs << ")");

    // Try to extract 1st element: axis min
    if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0), photon_energy_min)) {
        ERROR(errorPrefix << ", photon_energy_axis: 1st item must be a double (minimum photon energy)");
    }

    // Try to extract 2nd element: axis max
    if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1), photon_energy_max)) {
        ERROR(errorPrefix << ", photon_energy_axis: 2nd item must be a double (maximum photon energy)");
    }

    // Try to extract 3rd element: axis nbins
    if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2), photon_energy_nbins)) {
        ERROR(errorPrefix << ", photon_energy_axis: 3rd item must be an int (number of bins)");
    }

    // Check for  other keywords such as "logscale" and "edge_inclusive"
    photon_energy_logscale = false;
    photon_energy_edge_inclusive = false;
    for(unsigned int i=3; i<lenAxisArgs; i++) {
        std::string my_str("");
        PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
        if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
            photon_energy_logscale = true;
        else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
            photon_energy_edge_inclusive = true;
        else
            ERROR(errorPrefix << ": keyword `" << my_str << "` not understood");
    }

    // construct the list of photon_energies
    photon_energies.resize(photon_energy_nbins);
    delta_energies.resize(photon_energy_nbins);
    double emin=photon_energy_min, emax=photon_energy_max;
    if(photon_energy_logscale) {
        emin = log10(emin);
        emax = log10(emax);
    }
    spacing = (emax-emin)/photon_energy_nbins;
    for ( int i=0; i<photon_energy_nbins; i++) {
        photon_energies[i] = emin + (i+0.5)*spacing;
        if(photon_energy_logscale) {
            photon_energies[i] = pow(10.,photon_energies[i]);
            delta_energies[i]  = pow(10.,emin+i*spacing)*(pow(10.,spacing)-1.);
        } else {
            delta_energies[i]  = spacing;
        }
    }


    // here ends the reading of "photon energy axis"
    // ---------------------------------------------


    // Temporarily set the spatial min and max to the simulation box size
    spatial_min.resize( params.nDim_particle, 0. );
    spatial_max = params.grid_length;

    // get parameter "axes" that adds axes to the diagnostic
    // Each axis should contain several items:
    //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
    vector<PyObject*> pyAxes=PyTools::extract_pyVec("axes","DiagRadiationSpectrum",n_diag_rad_spectrum);

    // Create the Histogram object based on the extracted parameters above
    vector<string> excluded_axes(0);
    excluded_axes.push_back( "a" );
    excluded_axes.push_back( "b" );
    excluded_axes.push_back( "theta" );
    excluded_axes.push_back( "phi" );
    PyObject *deposited_quantity= PyString_FromString("dummy_radiation_spectrum");
    //histogram = HistogramFactory::create(params, NULL, pyAxes, species, patch, excluded_axes, errorPrefix);
    histogram = HistogramFactory::create(params, deposited_quantity, pyAxes, species, patch, excluded_axes, errorPrefix);

    // Get info on the spatial extent
    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
        if        ( histogram->axes[i]->type == "x" ) {
            spatial_min[0] = histogram->axes[i]->min;
            spatial_max[0] = histogram->axes[i]->max;
        } else if ( histogram->axes[i]->type == "y" ) {
            spatial_min[1] = histogram->axes[i]->min;
            spatial_max[1] = histogram->axes[i]->max;
        } else if ( histogram->axes[i]->type == "z" ) {
            spatial_min[2] = histogram->axes[i]->min;
            spatial_max[2] = histogram->axes[i]->max;
        }
    }

    // Calculate the size of the output array
    uint64_t total_size = photon_energy_nbins;
    for( unsigned int i=0; i<histogram->axes.size(); i++ )
        total_size *= histogram->axes[i]->nbins;
    if( total_size > 4294967296 ) // 2^32
        ERROR(errorPrefix << ": too many points (" << total_size << " > 2^32)");
    output_size = (unsigned int) total_size;

    // Output info on diagnostics
    if ( smpi->isMaster() ) {
        ostringstream mystream("");
        mystream.str("");
        mystream << species_names[0];
        for(unsigned int i=1; i<species_names.size(); i++)
            mystream << "," << species_names[i];
        MESSAGE(1,"Created RadiationSpectrum diagnostic #" << n_diag_rad_spectrum << ": species " << mystream.str());
        mystream.str("");
        mystream << "Photon energy from " << photon_energy_min << " to " << photon_energy_max << " in " << photon_energy_nbins << " steps";
        if( photon_energy_logscale       ) mystream << " [LOGSCALE] ";
        if( photon_energy_edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
        MESSAGE(2,mystream.str());
        for(unsigned int i=0; i<histogram->axes.size(); i++) {
            HistogramAxis * axis = histogram->axes[i];
            mystream.str("");
            mystream << "Axis " << axis->type << " from " << axis->min << " to " << axis->max << " in " << axis->nbins << " steps";
            if( axis->logscale       ) mystream << " [LOGSCALE] ";
            if( axis->edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
            MESSAGE(2,mystream.str());
        }

        // init HDF files (by master, only if it doesn't yet exist)
        mystream.str(""); // clear
        mystream << "RadiationSpectrum" << n_diag_rad_spectrum << ".h5";
        filename = mystream.str();
    }

} // END DiagnosticRadiationSpectrum::DiagnosticRadiationSpectrum


DiagnosticRadiationSpectrum::~DiagnosticRadiationSpectrum()
{
    delete timeSelection;
    delete flush_timeSelection;
} // END DiagnosticRadiationSpectrum::~DiagnosticRadiationSpectrum


// Called only by patch master of process master
void DiagnosticRadiationSpectrum::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if (!smpi->isMaster()) return;

    if( fileId_>0 ) return;

    if ( newfile ) {
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        // write all parameters as HDF5 attributes
        H5::attr(fileId_, "Version", string(__VERSION));
        H5::attr(fileId_, "time_average"  , time_average);
        // write all species
        ostringstream mystream("");
        mystream.str(""); // clear
        for (unsigned int i=0 ; i < species.size() ; i++)
            mystream << species[i] << " ";
        H5::attr(fileId_, "species", mystream.str());
        // write photon_energy_axis
        mystream.str(""); // clear
        mystream << "photon_energy_axis";
        string str1 = mystream.str();
        mystream.str(""); // clear
        mystream << photon_energy_min << " " << photon_energy_max << " "
        << photon_energy_nbins << " " << photon_energy_logscale << " " << photon_energy_edge_inclusive;
        string str2 = mystream.str();
        H5::attr(fileId_, str1, str2);
        // write each axis
        for (unsigned int iaxis=0 ; iaxis < histogram->axes.size() ; iaxis++) {
            mystream.str(""); // clear
            mystream << "axis" << iaxis;
            string str1 = mystream.str();
            mystream.str(""); // clear
            mystream << histogram->axes[iaxis]->type << " " << histogram->axes[iaxis]->min << " " << histogram->axes[iaxis]->max << " "
                     << histogram->axes[iaxis]->nbins << " " << histogram->axes[iaxis]->logscale << " " << histogram->axes[iaxis]->edge_inclusive << " [";
            for( unsigned int idim=0; idim<histogram->axes[iaxis]->coefficients.size(); idim++) {
                mystream << histogram->axes[iaxis]->coefficients[idim];
                if(idim<histogram->axes[iaxis]->coefficients.size()-1) mystream << ",";
            }
            mystream << "]";
            string str2 = mystream.str();
            H5::attr(fileId_, str1, str2);
        }
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    else {
        fileId_ = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
}


void DiagnosticRadiationSpectrum::closeFile()
{
    if (fileId_!=0) {
        H5Fclose(fileId_);
        fileId_ = 0;
    }

} // END closeFile


bool DiagnosticRadiationSpectrum::prepare( int timestep )
{
    // Get the previous timestep of the time selection
    int previousTime = timeSelection->previousTime(timestep);

    // Leave if the timestep is not the good one
    if (timestep - previousTime >= time_average) return false;

    // Allocate memory for the output array (already done if time-averaging)
    data_sum.resize(output_size);

    // if first time, erase output array
    if (timestep == previousTime)
        fill(data_sum.begin(), data_sum.end(), 0.);

    return true;

} // END prepare


// run one particle binning diagnostic
void DiagnosticRadiationSpectrum::run( Patch* patch, int timestep, SimWindow* simWindow )
{

    vector<int> int_buffer;
    vector<double> double_buffer;
    unsigned int npart, ndim = spatial_min.size();

    // Update spatial_min and spatial_max if needed
    for( unsigned int i=0; i<histogram->axes.size(); i++ ) {
        if( histogram->axes[i]->type == "moving_x" ) {
            spatial_min[0] = histogram->axes[i]->min + simWindow->getXmoved();
            spatial_max[0] = histogram->axes[i]->max + simWindow->getXmoved();
        }
    }

    // Verify that this patch is in a useful region for this diag
    for( unsigned int idim=0; idim<ndim; idim++ )
        if( patch->getDomainLocalMax(idim) < spatial_min[idim]
         || patch->getDomainLocalMin(idim) > spatial_max[idim] )
            return;

    // loop species
    for (unsigned int ispec=0 ; ispec < species.size() ; ispec++) {

        Species * s = patch->vecSpecies[species[ispec]];
        npart = s->particles->size();
        int_buffer   .resize(npart);
        double_buffer.resize(npart);

        fill(int_buffer.begin(), int_buffer.end(), 0);

        histogram->digitize  ( s, double_buffer, int_buffer, simWindow );

        // Sum the data into the data_sum
        // ------------------------------
        int ind;

        double gamma_inv, gamma, chi, xi, zeta, nu, cst;
        double two_third_ov_chi, increment0, increment;
        int iphoton_energy_max;

        double emin=photon_energy_min, emax=photon_energy_max;
        if (photon_energy_logscale) {
            emin = log10(emin);
            emax = log10(emax);
        }

        for (unsigned int ipart = 0 ; ipart < npart ; ipart++) {
            ind = int_buffer[ipart];
            if (ind<0) continue; // skip discarded particles
            ind *= photon_energy_nbins;

            // Get the quantum parameter
            chi = s->particles->chi(ipart);

            // Update the spectrum only if the quantum parameter is sufficiently high
            if (chi > 1e-5)
            {

                // Emitting particle energy (maximum of the spectrum)
                gamma = sqrt( 1. + pow(s->particles->momentum(0,ipart),2)
                                        + pow(s->particles->momentum(1,ipart),2)
                                        + pow(s->particles->momentum(2,ipart),2) );
                gamma_inv = 1./gamma;

                two_third_ov_chi = two_third/chi;

                increment0 = gamma_inv * s->particles->weight(ipart);

                // Comput the maximum iteration of the loop on bins
                if (photon_energy_logscale)
                {
                    iphoton_energy_max = int((log10(gamma) - emin)/spacing);
                }
                else
                {
                    iphoton_energy_max = int((gamma - emin)/spacing);
                }
                //iphoton_energy_max can not be greater than photon_energy_nbins
                iphoton_energy_max = min(iphoton_energy_max,photon_energy_nbins);


                // Loop on bins
                for (int i=0; i<iphoton_energy_max; i++) {

                    xi   = photon_energies[i] * gamma_inv;
                    //if ( xi<1.) {
                    zeta = xi/ (1.-xi);
                    nu   = two_third_ov_chi * zeta;
                    cst  = xi*zeta;
                    increment = increment0 * delta_energies[i]
                    //*           xi * ( RadiationTools::compute_f1_nu(nu) + cst*RadiationTools::compute_f2_nu(nu) );
                    *             xi * RadiationTools::compute_bessel_parts_radiated_power(nu,cst);
                    #pragma omp atomic
                    data_sum[ind+i] += increment;
                    //}
                }
            }
        }

    }

} // END run


// Now the data_sum has been filled
// if needed now, store result to hdf file
void DiagnosticRadiationSpectrum::write(int timestep, SmileiMPI* smpi)
{
    if ( !smpi->isMaster() ) return;

    if (timestep - timeSelection->previousTime() != time_average-1) return;

    double coeff;
    coeff = factor;
    // if time_average, then we need to divide by the number of timesteps
    if (time_average > 1) {
        coeff /= ((double)time_average);
    }
    for (unsigned int i=0; i<output_size; i++)
        data_sum[i] *= coeff;

    // make name of the array
    ostringstream mystream("");
    mystream.str("");
    mystream << "timestep" << setw(8) << setfill('0') << timestep;

    // write the array if it does not exist already
    if (! H5Lexists( fileId_, mystream.str().c_str(), H5P_DEFAULT ) ) {
        // Prepare array dimensions
        unsigned int naxes = histogram->axes.size() + 1;
        hsize_t dims[naxes];
        for( unsigned int iaxis=0; iaxis<naxes-1; iaxis++) dims[iaxis] = histogram->axes[iaxis]->nbins;
        dims[naxes-1] = photon_energy_nbins;
        // Create file space
        hid_t sid = H5Screate_simple(naxes, &dims[0], NULL);
        hid_t pid = H5Pcreate(H5P_DATASET_CREATE);
        // create dataset
        hid_t did = H5Dcreate(fileId_, mystream.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
        // write vector in dataset
        H5Dwrite(did, H5T_NATIVE_DOUBLE, sid, sid, H5P_DEFAULT, &data_sum[0]);
        // close all
        H5Dclose(did);
        H5Pclose(pid);
        H5Sclose(sid);
    }

    if( flush_timeSelection->theTimeIsNow(timestep) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

    // Clear the array
    clear();
    data_sum.resize(0);
} // END write


//! Clear the array
void DiagnosticRadiationSpectrum::clear() {
    data_sum.resize(0);
    vector<double>().swap( data_sum );
}


// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticRadiationSpectrum::getDiskFootPrint(int istart, int istop, Patch* patch)
{
    uint64_t footprint = 0;

    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore(istop) - timeSelection->howManyTimesBefore(istart);

    // Add necessary global headers approximately
    footprint += 1500;

    // Add necessary timestep headers approximately
    footprint += ndumps * 640;

    // Add size of each dump
    footprint += ndumps * (uint64_t)(output_size) * 8;

    return footprint;
}
