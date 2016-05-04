
#include <string>
#include <sstream>

#include "DiagnosticTrack.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int diagId, int speciesId ) :
nDim_particle(params.nDim_particle)
{
    diagId_ = diagId; // Warning, not the index of the DiagTrack, but of all local diags
    speciesId_ = speciesId;
    species = patch->vecSpecies[speciesId_];
    
    string nameSpec = "Ntot_"+ species->species_type;
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_INDEPENDENT);
    iter = 0;
    
    // Get the time selection from the particles
    timeSelection = species->particles->track_timeSelection;
    
    type_ = "Track";

}


// cloning constructor
DiagnosticTrack::DiagnosticTrack(DiagnosticTrack* track, Patch* patch)
{
    nDim_particle = track->nDim_particle;
    speciesId_    = track->speciesId_;
    species       = patch->vecSpecies[speciesId_];
    transfer      = track->transfer;
    iter          = track->iter;
    timeSelection = new TimeSelection(track->timeSelection);
    type_ = "Track";
}


DiagnosticTrack::~DiagnosticTrack()
{
}


void DiagnosticTrack::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    
    
    if ( newfile ) {
        // Create HDF5 file
        ostringstream hdf_filename("");
        hdf_filename << "TrackParticles_" << species->species_type  << ".h5" ;
        filename = hdf_filename.str();
        
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Define maximum size
        //hsize_t maxDimsPart[2] = {H5S_UNLIMITED, (hsize_t)nParticles};
        
        dims[0] = 0;//timeSelection->numberOfEvents(0, params.n_time);
        if ( !nbrParticles_ ) {
            ERROR("DiagTrack empty or number of Particles in diag is null");
        }
        else
            dims[1] = nbrParticles_;
        
        //hid_t file_space = H5Screate_simple(2, dims, NULL);
        // Define maximum size
        hsize_t maxDimsPart[2] = {H5S_UNLIMITED, (hsize_t)nbrParticles_};
        hid_t file_space = H5Screate_simple(2, dims, maxDimsPart);   
        
        // Create the overall dataset in the HDF5 file with a new chunk every timestep
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        hsize_t chunk_dims[2] = {1, (hsize_t)nbrParticles_};
        H5Pset_chunk(plist, 2, chunk_dims);
        
        // Create the datasets for x, y and z
        hid_t did;
        unsigned int nPosition = species->particles->Position.size();
        for (unsigned int i=0; i<nPosition; i++) {
            ostringstream namePos("");
            namePos << "Position-" << i;
            did = H5Dcreate(fileId_, namePos.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        
        // Create the datasets for px, py and pz
        unsigned int nMomentum = species->particles->Momentum.size();
        for (unsigned int i=0; i<nMomentum; i++) {
            ostringstream nameMom("");
            nameMom << "Momentum-" << i;
            did = H5Dcreate(fileId_, nameMom.str().c_str(), H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        // Create the datasets for weight, charge and ID
        did = H5Dcreate(fileId_, "Weight", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        did = H5Dcreate(fileId_, "Charge", H5T_NATIVE_SHORT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        did = H5Dcreate(fileId_, "Id", H5T_NATIVE_UINT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        
        H5Pclose(plist);
        H5Sclose(file_space);
        
        // Create the dataset for the time array
        hsize_t ntimes[1] = { 0 };
        hsize_t maxtimes[1] = { H5S_UNLIMITED };
        file_space = H5Screate_simple(1, ntimes, maxtimes);

        // Create the overall dataset in the HDF5 file with a new chunk every timestep
        plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        hsize_t tchunk_dims[1] = {1};
        H5Pset_chunk(plist, 1, tchunk_dims);

        did = H5Dcreate(fileId_, "Times", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        H5Pclose(plist);
        H5Sclose(file_space);
        
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        
    }
    else {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid);
        H5Pclose(pid);
   }

}


void DiagnosticTrack::closeFile()
{
    H5Fclose( fileId_ );
}


bool DiagnosticTrack::prepare( int timestep )
{
    if( ! timeSelection->theTimeIsNow(timestep) ) return false;
    
    dims[0] ++;
    
    // For all dataspace
    vector<string> datasets;
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        ostringstream namePos("");
        namePos << "Position-" << idim;
        datasets.push_back( namePos.str() );
    }
    for (int idim=0 ; idim<3 ; idim++) {
        ostringstream nameMom("");
        nameMom << "Momentum-" << idim;
        datasets.push_back( nameMom.str() );
    }
    datasets.push_back( "Weight" );
    datasets.push_back( "Charge" );
    datasets.push_back( "Id" );
    
    for (int idsets = 0 ; idsets<datasets.size() ; idsets++) {
        string name = datasets[idsets];
        hid_t did = H5Dopen( fileId_, name.c_str(), H5P_DEFAULT );
        // Increase the size of the array with the previously defined size
        H5Dset_extent(did, dims);
        H5Dclose(did);
    }
    
    hsize_t tdims[1] = { dims[0] };
    hid_t did = H5Dopen( fileId_, "Times", H5P_DEFAULT );
    H5Dset_extent(did, tdims);
    H5Dclose(did);

    // Write the current timestep
    vector<hsize_t> locator(1, iter);
    hsize_t onetime[1] = { 1 };
    hid_t mem_space = H5Screate_simple(1, onetime, NULL);
    append( fileId_, "Times", timestep, mem_space, 1, H5T_NATIVE_INT, locator );
    H5Sclose( mem_space );


    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    
    return true;
}


void DiagnosticTrack::run( Patch* patch, int timestep )
{
}


void DiagnosticTrack::write(int timestep)
{
    iter ++;
    
    int locNbrParticles = species->getNbrOfParticles();
    
    // Create the locator (gives locations where to store particles in the file)
    vector<hsize_t> locator (locNbrParticles*2);
    for(int i=0; i<locNbrParticles; i++) {
        locator[i*2  ] = iter-1;
        locator[i*2+1] = species->particles->id(i)-1; // because particles label Id starts at 1
    }
    
    // Specify the memory dataspace (the size of the local array)
    hsize_t count[2] = {1, (hsize_t)locNbrParticles};
    hid_t mem_space = H5Screate_simple(2, count, NULL);
    
    // For each dataspace (x, y, z, px, py, pz, weight, charge and ID), add the local
    // array to the HDF5 file
    ostringstream namePos("");
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        namePos.str("");
        namePos << "Position-" << idim;
        append( fileId_, namePos.str(), species->particles->position(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    ostringstream nameMom("");
    for (int idim=0 ; idim<3 ; idim++) {
        nameMom.str("");
        nameMom << "Momentum-" << idim;
        append( fileId_, nameMom.str(), species->particles->momentum(idim)[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    }
    append( fileId_, "Weight", species->particles->weight()[0], mem_space, locNbrParticles, H5T_NATIVE_DOUBLE, locator );
    append( fileId_, "Charge", species->particles->charge()[0], mem_space, locNbrParticles, H5T_NATIVE_SHORT , locator );
    append( fileId_, "Id"    , species->particles->id()    [0], mem_space, locNbrParticles, H5T_NATIVE_UINT  , locator );
    
    H5Sclose( mem_space );
    
    //MPI_Barrier(MPI_COMM_WORLD); // synchro to manage differently
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

}

template <class T>
void DiagnosticTrack::append( hid_t fid, string name, T & property,  hid_t  mem_space, int nParticles, hid_t type, vector<hsize_t> &locator) {
    
    // Open existing dataset
    hid_t did = H5Dopen( fid, name.c_str(), H5P_DEFAULT );
    
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


void DiagnosticTrack::setFileSplitting( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches )
{
    
    // 1 - Internal patches offset
    
    int localNbrParticles = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        DiagnosticTrack* diag = static_cast<DiagnosticTrack*>(vecPatches(ipatch)->localDiags[diagId_]);
        diag->species->particles->addIdOffsets(localNbrParticles);
        localNbrParticles += diag->species->getNbrOfParticles();
    }
    
    // 2 - MPI offset
    
    // Get the number of particles for each MPI
    int sz = smpi->getSize();
    std::vector<int> allNbrParticles(sz, 0);
    MPI_Allgather( &localNbrParticles, 1, MPI_INT, &allNbrParticles[0], 1, MPI_INT, MPI_COMM_WORLD );
    
    // Calculate the cumulative sum
    for (int irk=1 ; irk<sz ; irk++)
        allNbrParticles[irk] += allNbrParticles[irk-1];
    
    // Calculate the MPI offset
    int offset = 0;
    if( ! smpi->isMaster() ) offset = allNbrParticles[smpi->getRank()-1];
    
    // Apply the MPI offset
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        DiagnosticTrack* diag = static_cast<DiagnosticTrack*>(vecPatches(ipatch)->localDiags[diagId_]);
        diag->species->particles->addIdOffsets( offset );
        diag->setGlobalNbrParticles(allNbrParticles[sz-1]);
    }

}

