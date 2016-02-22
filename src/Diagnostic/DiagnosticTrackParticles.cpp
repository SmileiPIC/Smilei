#include "DiagnosticTrackParticles.h"

#include <iomanip>
#include <ostream>
#include <cmath>

#include "H5.h"
#include "Params.h"
#include "Patch.h"

using namespace std;

// constructor
DiagnosticTrackParticles::DiagnosticTrackParticles(Params& params, Patch* patch, Species* my_species) :
species(my_species),
nDim_particle(params.nDim_particle),
fid_(0)
{
    int locNbrParticles = species->getNbrOfParticles();
    hsize_t nParticles(0);
    int n = nParticles;

    string nameSpec="Ntot_"+species->species_type;

    // DRIVEN by PatchVector
    //nParticles = (hsize_t) diags->getScalar(nameSpec);

    //dims[0] = 0;
    //dims[1] = n;
    
    // if patch->isMaster()
    //     Create file

    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_INDEPENDENT);
    iter = 0;
    
    // Get the time selection from the particles
    timeSelection = species->particles->track_timeSelection;
}


// The master proc creates the HDF5 file and defines the dataspaces
//if ( patch->isMaster() ) {
void DiagnosticTrackParticles::createFile(int nParticles, Params &params) {
    // Create HDF5 file
    ostringstream filename("");
    filename << "TrackParticles_" << species->species_type  << ".h5" ;
    hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
    fid_ = H5Fcreate( filename.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        
//    // Write attribute: track_every
//    H5::attr(fid_, "every", species->particles->track_every);
        
    // Define maximum size
    hsize_t maxDimsPart[2] = {H5S_UNLIMITED, (hsize_t)nParticles};

//    dims[0] = params.n_time / species->particles->track_every + 1;
    dims[0] = timeSelection->numberOfEvents(0, params.n_time);
    dims[1] = nParticles;

    hid_t file_space = H5Screate_simple(2, dims, NULL);
        
    // Create the overall dataset in the HDF5 file with a new chunk every timestep
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    //H5Pset_layout(plist, H5D_CHUNKED);
    //H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
    hsize_t chunk_dims[2] = {1, (hsize_t)nParticles};
    //H5Pset_chunk(plist, 2, chunk_dims);
    // Create the datasets for x, y and z
    hid_t did;
    unsigned int nPosition = species->particles->Position.size();
    for (unsigned int i=0; i<nPosition; i++) {
	ostringstream namePos("");
	namePos << "Position-" << i;
	did = H5Dcreate(fid_, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
    }
        
    // Create the datasets for px, py and pz
    unsigned int nMomentum = species->particles->Momentum.size();
    for (unsigned int i=0; i<nMomentum; i++) {
	ostringstream nameMom("");
	nameMom << "Momentum-" << i;
	did = H5Dcreate(fid_, nameMom.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
	H5Dclose(did);
    }
        
    // Create the datasets for weight, charge and ID
    did = H5Dcreate(fid_, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dclose(did);
    did = H5Dcreate(fid_, "Charge", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dclose(did);
    did = H5Dcreate(fid_, "Id", H5T_NATIVE_UINT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dclose(did);
        
    H5Pclose(plist);
    H5Sclose(file_space);
    
    // Create the dataset for the time array
    plist = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t ntimes[1] = { dims[0] };
    file_space = H5Screate_simple(1, ntimes, NULL);
    did = H5Dcreate(fid_, "Times", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dclose(did);
    H5Pclose(plist);
    H5Sclose(file_space);
    

    H5Fclose( fid_ );
}

void DiagnosticTrackParticles::open() {
    // Define the HDF5 MPI file access
    file_access_ = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio( file_access_, MPI_COMM_WORLD, MPI_INFO_NULL );
        
    // Open the HDF5 file
    ostringstream filename("");
    filename << "TrackParticles_" << species->species_type  << ".h5" ;
    fid_ = H5Fopen( filename.str().c_str(), H5F_ACC_RDWR, file_access_);
}

void DiagnosticTrackParticles::setFile( hid_t file_access, hid_t fid ) {
    file_access_ = file_access;
    fid_ = fid;
}

void DiagnosticTrackParticles::close() {
    H5Pclose( file_access_ );
    H5Fclose( fid_ );
}

void DiagnosticTrackParticles::run(int time) {
    
//    if ( time % species->particles->track_every != 0) return;
    if( ! timeSelection->theTimeIsNow(time) ) return;
    
    iter ++;

    int locNbrParticles = species->getNbrOfParticles();
    
    // Create the locator (gives locations where to store particles in the file)
    vector<hsize_t> locator (locNbrParticles*2);
    for(int i=0; i<locNbrParticles; i++) {
        locator[i*2  ] = iter-1;
        locator[i*2+1] = species->particles->id(i)-1; // because particles label Id starts at 1
    }
    // Now we increase the array size for the new timestep (new chunk)
    // It is not applied to the HDF5 file yet
    dims[0] ++;

    
    // Specify the memory dataspace (the size of the local array)
    hsize_t count[2] = {1, (hsize_t)locNbrParticles};
    hid_t mem_space = H5Screate_simple(2, count, NULL);

    // For each dataspace (x, y, z, px, py, pz, weight, charge and ID), add the local
    // array to the HDF5 file -> see function appendTestParticles() in SmileiIO.h
    ostringstream namePos("");
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        namePos.str("");
        namePos << "Position-" << idim;
        append( fid_, namePos.str(), species->particles->position(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    ostringstream nameMom("");
    for (int idim=0 ; idim<3 ; idim++) {
        nameMom.str("");
        nameMom << "Momentum-" << idim;
        append( fid_, nameMom.str(), species->particles->momentum(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    append( fid_, "Weight", species->particles->weight()[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    append( fid_, "Charge", species->particles->charge()[0], mem_space, locNbrParticles, H5T_NATIVE_SHORT , locator );
    append( fid_, "Id"    , species->particles->id()    [0], mem_space, locNbrParticles, H5T_NATIVE_UINT  , locator );
    
    H5Sclose( mem_space );
    
    // Write the current timestep
    locator.resize(1, iter-1);
    hsize_t onetime[1] = { 1 };
    mem_space = H5Screate_simple(1, onetime, NULL);
    append( fid_, "Times", time, mem_space, 1, H5T_NATIVE_INT, locator );
    
    //MPI_Barrier(MPI_COMM_WORLD); // synchro to manage differently
    H5Fflush( fid_, H5F_SCOPE_GLOBAL );


}


template <class T>
void DiagnosticTrackParticles::append( hid_t fid, string name, T & property,  hid_t  mem_space, int nParticles, hid_t type, vector<hsize_t> &locator) {
    
    // Open existing dataset
    hid_t did = H5Dopen( fid, name.c_str(), H5P_DEFAULT );
    // Increase the size of the array with the previously defined size

    //H5Dset_extent(did, dims);

    // Get the extended file space
    hid_t file_space = H5Dget_space(did);
    
    // Select locations that this proc will write
    if(nParticles>0)
        H5Sselect_elements( file_space, H5S_SELECT_SET, nParticles, &locator[0] );
    else
        H5Sselect_none(file_space);
    
    // Write
    H5Dwrite( did, type, mem_space , file_space , transfer, &property );
    
    H5Sclose(file_space);
    H5Dclose(did);
}

