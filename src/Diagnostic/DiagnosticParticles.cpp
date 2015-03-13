#include "DiagnosticParticles.h"

#include <iomanip>
#include <ostream>
#include <cmath>

#include "PicParams.h"
#include "DiagParams.h"

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
    for (int iaxis=0 ; iaxis < axes.size() ; iaxis++)
        output_size *= axes[iaxis]->nbins;
    
    // if necessary, resize the output array
    if (time_average>1)
        data_sum.resize(output_size);
    
    MESSAGE("Created particle diagnostic #" << ID);
    
}

// destructor
DiagnosticParticles::~DiagnosticParticles()
{
}


// -------------------
// Some HDF5 overlays
// -------------------
// write a string as an attribute
void H5_attr_string(int fileId, string attribute_name, string attribute_value) {
    hid_t sid  = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, attribute_value.size());
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    hid_t aid = H5Acreate(fileId, attribute_name.c_str(), atype, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, atype, attribute_value.c_str());
    H5Aclose(aid);
    H5Sclose(sid);
    H5Tclose(atype);
}
// write an unsigned int as an attribute
void H5_attr_uint(int fileId, string attribute_name, unsigned int attribute_value) {
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(fileId, attribute_name.c_str(), H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_UINT, &attribute_value);
    H5Sclose(sid);
    H5Aclose(aid);
}
// write an int as an attribute
void H5_attr_int(int fileId, string attribute_name, int attribute_value) {
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(fileId, attribute_name.c_str(), H5T_NATIVE_INT, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_INT, &attribute_value);
    H5Sclose(sid);
    H5Aclose(aid);
}
// write a double as an attribute
void H5_attr_double(int fileId, string attribute_name, double attribute_value) {
    hid_t sid = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate(fileId, attribute_name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_DOUBLE, &attribute_value);
    H5Sclose(sid);
    H5Aclose(aid);
}

// -------------------




// Declare static variables here
vector<DiagnosticParticles*> DiagnosticParticles::vecDiagnosticParticles;

// close all the files
void DiagnosticParticles::closeAll()
{
    
    int n = vecDiagnosticParticles.size();
    for (int i=0; i<n; i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->close();
    
}

// close the hdf file
void DiagnosticParticles::close()
{
    
    if (fileId != 0) H5Fclose(fileId);
    
}

// run all the particle diagnostics
void DiagnosticParticles::runAll(int timestep, vector<Species*> &vecSpecies, SmileiMPI* smpi)
{
    
    int n = vecDiagnosticParticles.size();
    for (int i=0; i<n; i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->run(timestep, vecSpecies, smpi);
    
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
    
    // if first timestep, init HFD files
    if (timestep ==0) {
        // the master creates the hdf file
        if (smpi->isMaster() ) {
            mystream.str("");
            mystream << "ParticleDiagnostic" << diagnostic_id << ".h5";
            fileId = H5Fcreate( mystream.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
            // write version
            string ver(__VERSION);
            H5_attr_string(fileId, "Version", ver);
            // write output
            H5_attr_string(fileId, "output" , output);
            // write every
            H5_attr_uint  (fileId, "every"  , every);
            // write time_average
            H5_attr_uint  (fileId, "time_average"  , time_average);
            // write all species
            mystream.str(""); // clear
            for (int i=0 ; i < species.size() ; i++)
                mystream << species[i] << " ";
            H5_attr_string(fileId, "species", mystream.str().c_str());
            // write each axis
            for (int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
                mystream.str(""); // clear
                mystream << "axis" << iaxis;
                str1 = mystream.str();
                mystream.str(""); // clear
                mystream << axes[iaxis]->type << " " << axes[iaxis]->min << " " << axes[iaxis]->max << " "
                         << axes[iaxis]->nbins << " " << axes[iaxis]->logscale << " " << axes[iaxis]->edge_inclusive;
                str2 = mystream.str();
                H5_attr_string(fileId, str1, str2);
            }
        }
    }
    
    // skip the routine if the timestep is not the good one
    if (timestep % every >= time_average) return;
    
    // create output array if not already available
    if (time_average == 1) { // if no time-average, then we need an array.
        data_sum.resize(output_size);
    }
    
    // if first time, erase output array
    if (timestep % every == 0)
        fill(data_sum.begin(), data_sum.end(), 0.);
    
    // loop species
    for (int ispec=0 ; ispec < species.size() ; ispec++) {
        
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
        for (int ibin=0 ; ibin<nbins ; ibin++) {
            
            bmin = s->bmin[ibin];
            bmax = s->bmax[ibin];
            
            // 1 - loop on the different axes requested and compute the output index of each particle
            // --------------------------------------------------------------------------------------
            for (int iaxis=0 ; iaxis < axes.size() ; iaxis++) {
                
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
            
            int imin=10000000, imax=0;
            for (int i=0; i<index_array.size(); i++) {
                if (index_array[i]<imin) imin = index_array[i];
                if (index_array[i]>imax) imax = index_array[i];
            }
            
            // 2 - prepare the data to output
            // ------------------------------
            if      (output == "density")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart];
            
            else if (output == "current_density_x")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (*px)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "current_density_y")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (*py)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
            else if (output == "current_density_z")
                for (int ipart = bmin ; ipart < bmax ; ipart++)
                    data_array[ipart] = (*w)[ipart] * (*pz)[ipart] / sqrt( 1. + pow((*px)[ipart],2) + pow((*py)[ipart],2) + pow((*pz)[ipart],2) );
            
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
                coeff = 1./((double)output_size);
                for (int i=0; i<output_size; i++)
                    data_sum[i] *= coeff;
            }
            
            // create dataspace for 1D array with "output_size" number of elements
            hsize_t dims = output_size;
            hid_t sid = H5Screate_simple(1, &dims, NULL);
            hid_t pid = H5Pcreate(H5P_DATASET_CREATE); // property list
            // make name of the dataset
            mystream.str("");
            mystream << "timestep" << setw(8) << setfill('0') << timestep;
            // create dataset 
            hid_t did = H5Dcreate(fileId, mystream.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, pid, H5P_DEFAULT);
            // write data_sum in dataset
            H5Dwrite(did, H5T_NATIVE_DOUBLE, sid, sid, H5P_DEFAULT, &data_sum[0]);
            // close all
            H5Dclose(did);
            H5Pclose(pid);
            H5Sclose(sid);
        }
        
    }
    
    // delete temporary stuff
    if (time_average == 1) data_sum.resize(0);
    
}


