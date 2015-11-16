#include "DiagnosticTrackParticles.h"

#include <iomanip>
#include <ostream>
#include <cmath>

#include "H5.h"
#include "Params.h"

using namespace std;

// constructor
DiagnosticTrackParticles::DiagnosticTrackParticles(Params& params, SmileiMPI* smpi, Species* my_species) :
species(my_species),
nDim_particle(params.nDim_particle)
{

    int locNbrParticles = species->getNbrOfParticles();
    hsize_t nParticles = (hsize_t) smpi->globalNbrParticles(species, locNbrParticles);
    
    // Define the initial array dimensions
    int n = nParticles;
    smpi->bcast(n); // necessary to broadcast nParticles because only master has it
    dims[0] = 0;
    dims[1] = n;
    
    // The master proc creates the HDF5 file and defines the dataspaces
    if ( smpi->isMaster() ) {
        
        // Create HDF5 file
        ostringstream filename("");
        filename << params.output_dir << "/TrackParticles_" << species->species_type  << ".h5" ;
        hid_t fid = H5Fcreate( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        
        // Write attribute: dump_every
        H5::attr(fid, "every", species->particles.dump_every);
        
        // Define maximum size
        hsize_t maxDimsPart[2] = {H5S_UNLIMITED, nParticles};
        hid_t file_space = H5Screate_simple(2, dims, maxDimsPart);
        
        // Create the overall dataset in the HDF5 file with a new chunk every timestep
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        hsize_t chunk_dims[2] = {1, nParticles};
        H5Pset_chunk(plist, 2, chunk_dims);
        
        // Create the datasets for x, y and z
        hid_t did;
        unsigned int nPosition = species->particles.Position.size();
        for (unsigned int i=0; i<nPosition; i++) {
            ostringstream namePos("");
            namePos << "Position-" << i;
            did = H5Dcreate(fid, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        // Create the datasets for px, py and pz
        unsigned int nMomentum = species->particles.Momentum.size();
        for (unsigned int i=0; i<nMomentum; i++) {
            ostringstream nameMom("");
            nameMom << "Momentum-" << i;
            did = H5Dcreate(fid, nameMom.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        // Create the datasets for weight, charge and ID
        did = H5Dcreate(fid, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        did = H5Dcreate(fid, "Charge", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        did = H5Dcreate(fid, "Id", H5T_NATIVE_UINT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        
        H5Pclose(plist);
        H5Sclose(file_space);
        H5Fclose( fid );
    }
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE);

}


void DiagnosticTrackParticles::run(int time, SmileiMPI* smpi) {
    
    if ( time % species->particles.dump_every != 0) return;

    int locNbrParticles = species->getNbrOfParticles();
    
    // We increase the array size for the new timestep (new chunk)
    // It is not applied to the HDF5 file yet
    dims[0] ++;
    
    // Create the locator (gives locations where to store particles in the file)
    vector<hsize_t> locator (locNbrParticles*2);
    for(int i=0; i<locNbrParticles; i++) {
        locator[i*2  ] = dims[0]-1;
        locator[i*2+1] = species->particles.id(i)-1;
    }
    
    // Define the HDF5 MPI file access
    hid_t file_access = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio( file_access, MPI_COMM_WORLD, MPI_INFO_NULL );
    
    // Specify the memory dataspace (the size of the local array)
    hsize_t count[2] = {1, (hsize_t)locNbrParticles};
    mem_space = H5Screate_simple(2, count, NULL);
    
    // Open the HDF5 file
    ostringstream filename("");
    filename << "TrackParticles_" << species->species_type  << ".h5" ;
    hid_t fid = H5Fopen( filename.str().c_str(), H5F_ACC_RDWR, file_access);

    // For each dataspace (x, y, z, px, py, pz, weight, charge and ID), add the local
    // array to the HDF5 file -> see function appendTestParticles() in SmileiIO.h
    ostringstream namePos("");
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        namePos.str("");
        namePos << "Position-" << idim;
        append( fid, namePos.str(), species->particles.position(idim), locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    ostringstream nameMom("");
    for (int idim=0 ; idim<3 ; idim++) {
        nameMom.str("");
        nameMom << "Momentum-" << idim;
        append( fid, nameMom.str(), species->particles.momentum(idim), locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    append( fid, "Weight", species->particles.weight(), locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    append( fid, "Charge", species->particles.charge(), locNbrParticles, H5T_NATIVE_SHORT , locator );
    append( fid, "Id"    , species->particles.id()    , locNbrParticles, H5T_NATIVE_UINT  , locator );
    
    smpi->barrier(); // sometimes needed because all procs must not close the file too early
    H5Fflush( fid, H5F_SCOPE_GLOBAL );
    H5Pclose( file_access );
    H5Fclose( fid );
    
}


template <class T>
void DiagnosticTrackParticles::append( hid_t fid, string name, std::vector<T> property, int nParticles, hid_t type, vector<hsize_t> &locator) {
    
    // Open existing dataset
    hid_t did = H5Dopen( fid, name.c_str(), H5P_DEFAULT );
    // Increase the size of the array with the previously defined size
    H5Dset_extent(did, dims);
    // Get the extended file space
    hid_t file_space = H5Dget_space(did);
    
    // Select locations that this proc will write
    if(nParticles>0)
        H5Sselect_elements( file_space, H5S_SELECT_SET, nParticles, &locator[0] );
    else
        H5Sselect_none(file_space);
    
    // Write
    H5Dwrite( did, type, mem_space , file_space , transfer, &(property[0]) );
    
    H5Sclose(file_space);
    H5Dclose(did);
}

