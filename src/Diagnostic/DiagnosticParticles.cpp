#include "DiagnosticParticles.h"

#include <iomanip>
#include <ostream>
#include <cmath>

#include "Params.h"
#include "H5.h"

using namespace std;

// constructor
DiagnosticParticles::DiagnosticParticles(unsigned int ID, string output_, unsigned int every_, unsigned int time_average_, vector<unsigned int> species_, vector<DiagnosticParticlesAxis*> axes_) :
    fileId(0)
{

    diagnostic_id = ID;
    output = output_;
    every = every_;
    time_average = time_average_;
    species = species_;
    axes = axes_;
    
    // calculate the total size of the output array
    output_size = 1;
    for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++)
        output_size *= axes[iaxis]->nbins;
    
    // if necessary, resize the output array
    if (time_average>1)
        data_sum.resize(output_size);
    
    // Output info on diagnostics
    ostringstream mystream("");
    mystream.str("");
    mystream << species[0];
    for(unsigned int i=1; i<species.size(); i++)
        mystream << "," << species[i];
    MESSAGE(2,"Created particle diagnostic #" << ID << ": species " << mystream.str());
    DiagnosticParticlesAxis *a;
    for(unsigned int i=0; i<axes.size(); i++) {
        a = axes[i];
        mystream.str("");
        mystream << "\t\t\tAxis " << a->type << " from " << a->min << " to " << a->max << " in " << a->nbins << " steps";
        if( a->logscale       ) mystream << " [LOGSCALE] ";
        if( a->edge_inclusive ) mystream << " [EDGE INCLUSIVE]";
        MESSAGE(mystream.str());
    }
    
    
}

// destructor
DiagnosticParticles::~DiagnosticParticles()
{
}

// close the hdf file
void DiagnosticParticles::close()
{
    
    if (fileId != 0) H5Fclose(fileId);
    
}

// run one particle diagnostic
void DiagnosticParticles::run(int timestep, vector<Species*>& vecSpecies, SmileiMPI* smpi)
{
    
    Species *s;
    Particles *p;
    vector<int> index_array;
    vector<double> *x, *y, *z, *px, *py, *pz, *w, axis_array, data_array;
    vector<short> *q;
    int nbins = vecSpecies[0]->bmin.size(); // number of bins in the particles binning (openMP)
    int bmin, bmax, axissize, ind;
    double axismin, axismax, mass, coeff;
    string axistype, str1, str2;
    ostringstream mystream("");
    
    // init HFD files (by master, only if it doesn't yet exist)
    if (smpi->isMaster() && !fileId) {
        mystream.str("");
        mystream << "ParticleDiagnostic" << diagnostic_id << ".h5";
        fileId = H5Fcreate( mystream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        // write all parameters as HDF5 attributes
        string ver(__VERSION);
        H5::attr(fileId, "Version", ver);
        H5::attr(fileId, "output" , output);
        H5::attr(fileId, "every"  , every);
        H5::attr(fileId, "time_average"  , time_average);
        // write all species
        mystream.str(""); // clear
        for (unsigned int i=0 ; i < species.size() ; i++)
            mystream << species[i] << " ";
        H5::attr(fileId, "species", mystream.str().c_str());
        // write each axis
        for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
            mystream.str(""); // clear
            mystream << "axis" << iaxis;
            str1 = mystream.str();
            mystream.str(""); // clear
            mystream << axes[iaxis]->type << " " << axes[iaxis]->min << " " << axes[iaxis]->max << " "
            << axes[iaxis]->nbins << " " << axes[iaxis]->logscale << " " << axes[iaxis]->edge_inclusive;
            str2 = mystream.str();
            H5::attr(fileId, str1, str2);
        }
    }
    
    // skip the routine if the timestep is not the good one
    if (timestep % every >= time_average) return;
    
    // Allocate memory for the output array (already done if time-averaging)
    if (time_average <= 1)
        data_sum.resize(output_size);
    
    // if first time, erase output array
    if (timestep % every == 0)
        fill(data_sum.begin(), data_sum.end(), 0.);
    
    // loop species
    for (unsigned int ispec=0 ; ispec < species.size() ; ispec++) {
        
        // make shortcuts
        s = vecSpecies[species[ispec]]; // current species
        p = &(s->particles);            // current particles array
        x  = &(p->Position[0]);         // -+
        y  = &(p->Position[1]);         //  |- position
        z  = &(p->Position[2]);         // -+
        px = &(p->Momentum[0]);         // -+
        py = &(p->Momentum[1]);         //  |- momentum
        pz = &(p->Momentum[2]);         // -+
        q  = &(p->Charge);              // charge
        w  = &(p->Weight);              // weight
        mass = s->species_param.mass;   // mass
        
        axis_array .resize(p->size()); // array to store particle axis data
        index_array.resize(p->size()); // array to store particle output index
        data_array .resize(p->size()); // array to store particle output data
        
        fill(index_array.begin(), index_array.end(), 0);
        
        // loop each openMP bin
        //! \todo Make OpenMP parallelization
        for (int ibin=0 ; ibin<nbins ; ibin++) {
            
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            
            // 1 - loop on the different axes requested and compute the output index of each particle
            // --------------------------------------------------------------------------------------
            for (unsigned int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
                
                axismin  = axes[iaxis]->min;
                axismax  = axes[iaxis]->max;
                axissize = axes[iaxis]->nbins;
                axistype = axes[iaxis]->type;
                
                // first loop on particles to store the indexing (axis) quantity
                
                if      (axistype == "x"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*x)[ipart];
                
                else if (axistype == "y"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*y)[ipart];
                
                else if (axistype == "y"     )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*z)[ipart];
                
                else if (axistype == "px"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*px)[ipart];
                
                else if (axistype == "py"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*py)[ipart];
                
                else if (axistype == "pz"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (*pz)[ipart];
                
                else if (axistype == "p"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * sqrt(pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2));
                
                else if (axistype == "gamma" )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "ekin" )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = mass * (sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) ) - 1.);
                
                else if (axistype == "vx"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*px)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "vy"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*py)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "vz"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (*pz)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                        
                else if (axistype == "v"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = pow( 1. + 1./(pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2)) , -0.5);
                
                else if (axistype == "vperp2"    )
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (pow((*py)[ipart],2) + pow((*pz)[ipart],2)) / (1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
                
                else if (axistype == "charge")
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = (double) (*q)[ipart];
                
                else    ERROR("In particle diagnostics, axis `" << axistype << "` unknown");
                
                // Now, axis_array points to the array that contains the particles data (for indexing)
                
                // if log scale
                if (axes[iaxis]->logscale) {
                    // then loop again and convert to log
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        axis_array[ipart] = log10(axis_array[ipart]);
                    // also convert other quantities
                    axismin = log10(axismin);
                    axismax = log10(axismax);
                }
                
                // The indexes are "reshaped" in one dimension.
                // For instance, in 3d, the index has the form  i = i3 + n3*( i2 + n2*i1 )
                // Here we do the multiplication by n3 or n2 (etc.)
                if (iaxis>0) {
                    for (int ipart = bmin ; ipart < bmax ; ipart++)
                        index_array[ipart] *= axissize;
                }
                
                // loop again on the particles and calculate index
                // This is separated in two cases: edge_inclusive and edge_exclusive
                coeff = ((double)axissize)/(axismax-axismin);
                if (!axes[iaxis]->edge_inclusive) { // if the particles out of the "box" must be excluded
                    
                    for (int ipart = bmin ; ipart < bmax ; ipart++) {
                        // skip already discarded particles
                        if (index_array[ipart] < 0) continue; 
                        
                        // calculate index
                        ind = floor( (axis_array[ipart]-axismin) * coeff );
                        
                        // index valid only if in the "box"
                        if (ind >= 0  &&  ind < axissize) index_array[ipart] += ind;
                        else index_array[ipart] = -1; // discard particle
                    }
                    
                } else { // if the particles out of the "box" must be included
                    
                    for (int ipart = bmin ; ipart < bmax ; ipart++) {
                        // skip already discarded particles
                        if (index_array[ipart] < 0) continue; 
                        
                        // calculate index
                        ind = floor( (axis_array[ipart]-axismin) * coeff );
                        
                        // move out-of-range indexes back into range
                        if (ind < 0) ind = 0;
                        if (ind >= axissize) ind = axissize-1;
                        index_array[ipart] += ind;
                    }
                    
                }
                
            } // loop axes
            
            // 2 - prepare the data to output
            // ------------------------------
            if      (output == "density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart];
            
            if      (output == "charge_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]);
            
            else if (output == "current_density_x")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*px)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "current_density_y")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*py)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "current_density_z")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (double)((*q)[ipart]) * (*pz)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "p_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * sqrt(pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2));
            
            else if (output == "px_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*px)[ipart];
            
            else if (output == "py_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*py)[ipart];
            
            else if (output == "pz_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * (*pz)[ipart];
            
            else if (output == "pxvx_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * pow((*px)[ipart],2)/ sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "pyvy_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * pow((*py)[ipart],2)/ sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "pzvz_density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = mass * (*w)[ipart] * pow((*pz)[ipart],2)/ sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            // 3 - sum the data into the data_sum according to the indexes
            // ---------------------------------------------------------------
            for (int ipart = bmin ; ipart < bmax ; ipart++) {
                ind = index_array[ipart];
                if (ind<0) continue; // skip discarded particles
                data_sum[ind] += data_array[ipart];
            }
            
        } // loop openMP bins
        
    } // loop species
    
    // Now the data_sum has been filled
    
    // if needed now, store result to hdf file
    if (timestep % every == time_average-1) {
        
        // sum the outputs from each MPI partition
        MPI_Reduce( smpi->isMaster()?MPI_IN_PLACE:&data_sum[0], &data_sum[0], output_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        
        if (fileId > 0) { // only the master has fileId>0
            // if time_average, then we need to divide by the number of timesteps
            if (time_average > 1) {
                coeff = 1./((double)time_average);
                for (int i=0; i<output_size; i++)
                    data_sum[i] *= coeff;
            }
            // make name of the array
            mystream.str("");
            mystream << "timestep" << setw(8) << setfill('0') << timestep;
            // write the array
            H5::vect(fileId, mystream.str(), data_sum);
            H5Fflush(fileId, H5F_SCOPE_GLOBAL);
        }
        
    }
    
    // delete temporary stuff
    if (time_average == 1) data_sum.resize(0);
    
}


