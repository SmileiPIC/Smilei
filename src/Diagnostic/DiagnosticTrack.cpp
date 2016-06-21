
#include <string>
#include <sstream>

#include "DiagnosticTrack.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int diagId, int speciesId ) :
nDim_particle(params.nDim_particle),
IDs_done( params.restart )
{
    diagId_ = diagId; // Warning, not the index of the DiagTrack, but of all local diags
    speciesId_ = speciesId;
    Species* species = patch->vecSpecies[speciesId_];
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE);
    
    // Get the time selection from the particles
    timeSelection = species->particles->track_timeSelection;
    
    // Create the filename
    ostringstream hdf_filename("");
    hdf_filename << "TrackParticles_" << species->species_type  << ".h5" ;
    filename = hdf_filename.str();
    
    // Create a list of the necessary datasets
    datasets.push_back( "Id" );
    datatypes.push_back( H5T_NATIVE_UINT );
    datasets.push_back( "Charge" );
    datatypes.push_back( H5T_NATIVE_SHORT );
    datasets.push_back( "Weight" );
    datatypes.push_back( H5T_NATIVE_DOUBLE );
    datasets.push_back( "Momentum-0" );
    datatypes.push_back( H5T_NATIVE_DOUBLE );
    datasets.push_back( "Momentum-1" );
    datatypes.push_back( H5T_NATIVE_DOUBLE );
    datasets.push_back( "Momentum-2" );
    datatypes.push_back( H5T_NATIVE_DOUBLE );
    ostringstream name;
    for (int idim=0 ; idim<nDim_particle ; idim++) {
        name.str("");
        name << "Position-" << idim;
        datasets.push_back( name.str() );
        datatypes.push_back( H5T_NATIVE_DOUBLE );
    }
    
    type_ = "Track";
}

DiagnosticTrack::~DiagnosticTrack()
{
}


void DiagnosticTrack::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    
    if ( newfile ) {
        // Create HDF5 file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Define maximum size
        //hsize_t maxDimsPart[2] = {H5S_UNLIMITED, (hsize_t)nParticles};
        
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
        
        // Create the datasets
        for (unsigned int i=0; i<datasets.size(); i++) {
            hid_t did = H5Dcreate(fileId_, datasets[i].c_str(), datatypes[i], file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            H5Dclose(did);
        }
        
        H5Pclose(plist);
        H5Sclose(file_space);
        
        // Create the dataset for the time array
        hsize_t ntimes[1] = { 0 };
        hsize_t maxtimes[1] = { H5S_UNLIMITED };
        file_space = H5Screate_simple(1, ntimes, maxtimes);
        plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        hsize_t chunks[1] = {1};
        H5Pset_chunk(plist, 1, chunks);
        hid_t did = H5Dcreate(fileId_, "Times", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        H5Pclose(plist);
        H5Sclose(file_space);
        
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        
    }
    else {
        // Open the file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid);
        
        // Open the "Id" group to get the current dimensions
        hid_t did = H5Dopen( fileId_, "Id", H5P_DEFAULT );
        hid_t sid = H5Dget_space( did );
        H5Sget_simple_extent_dims(sid, &dims[0], NULL );
        
        // Close all things but not the file
        H5Sclose(sid);
        H5Dclose(did);
        H5Pclose(pid);
   }

}


void DiagnosticTrack::closeFile()
{
    if(fileId_>0) {
        H5Fclose( fileId_ );
        fileId_=0;
    }
}


void DiagnosticTrack::init(SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // Set the IDs of the particles
    if( ! IDs_done ) {
        
        // 1 - Internal patches offset
        
        int localNbrParticles = 0;
        for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            vecPatches(ipatch)->vecSpecies[speciesId_]->particles->addIdOffsets(localNbrParticles);
            localNbrParticles += vecPatches(ipatch)->vecSpecies[speciesId_]->getNbrOfParticles();
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
            vecPatches(ipatch)->vecSpecies[speciesId_]->particles->addIdOffsets( offset );
        }
        nbrParticles_ = allNbrParticles[sz-1];
        
        IDs_done = true;
    }
    
    // define the initial arrays sizes
    dims[0] = 0;
    if ( !nbrParticles_ )
        ERROR("DiagTrack empty or number of Particles in diag is null");
    dims[1] = nbrParticles_;
    
}


bool DiagnosticTrack::prepare( int timestep )
{
    if( ! timeSelection->theTimeIsNow(timestep) ) return false;
    
    return true;

}


void DiagnosticTrack::run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep )
{
    
    // Add a new timestep to the dimension of the arrays
    dims[0]++;
    
    // Obtain the particle partition of all the patches in this MPI
    patch_start.resize( vecPatches.size() );
    int nParticles = 0;
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        patch_start[ipatch] = nParticles;
        nParticles += vecPatches(ipatch)->vecSpecies[speciesId_]->getNbrOfParticles();
    }
    
    // Build the "locator", an array indicating where each particle goes in the final array
    vector<hsize_t> locator (nParticles*2);
    for (int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        Particles* particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
        int np=particles->size(), i=0, j=patch_start[ipatch];
        while( i<np ) {
            locator[j*2  ] = dims[0]-1;
            locator[j*2+1] = particles->id(i)-1; // because particles label Id starts at 1
            i++; j++;
        }
    }
    
    // Specify the memory dataspace (the size of the local array)
    hsize_t count[1] = {(hsize_t)nParticles};
    mem_space = H5Screate_simple(1, count, NULL);
    
    // For each dataset
    for( int idset=0; idset<datasets.size(); idset++) {
        
        // Fill the buffer for the current patch
        unsigned int patch_nParticles, i, j;
        Particles* particles;
        if( idset == 0 ) { // Id
            data_uint.resize( nParticles );
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_uint[j] = particles->id(i);
                    i++; j++;
                }
            }
        } else if( idset == 1 ) { // Charge
            data_short.resize( nParticles );
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_short[j] = particles->charge(i);
                    i++; j++;
                }
            }
        } else if( idset == 2 ) { // Weight
            data_double.resize( nParticles );
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_double[j] = particles->weight(i);
                    i++; j++;
                }
            }
        } else if( idset < 6 ) { // Momentum
            data_double.resize( nParticles );
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_double[j] = particles->momentum(idset-3, i);
                    i++; j++;
                }
            }
        } else { // Position
            data_double.resize( nParticles );
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_double[j] = particles->position(idset-6, i);
                    i++; j++;
                }
            }
        }
        
        // Open existing dataset
        hid_t did = H5Dopen( fileId_, datasets[idset].c_str(), H5P_DEFAULT );
        // Update the size of this dataset for the new timestep
        H5Dset_extent(did, dims);
        
        // Get the extended file space
        hid_t file_space = H5Dget_space(did);
        
        // Select locations that this proc will write
        if(nParticles>0)
            H5Sselect_elements( file_space, H5S_SELECT_SET, nParticles, &locator[0] );
        else
            H5Sselect_none(file_space);
        
        // Write
        if( idset == 0 ) {
            H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_uint[0] );
        } else if( idset == 1 ) {
            H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_short[0] );
        } else {
            H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_double[0] );
        }
        
        H5Sclose(file_space);
        H5Dclose(did);
        
    }
    H5Sclose( mem_space );
    
    // Update the size of the times dataset for the new timestep
    hsize_t tdims[1] = { dims[0] };
    hid_t did = H5Dopen( fileId_, "Times", H5P_DEFAULT );
    H5Dset_extent(did, tdims);
    // Select only the last element of the array
    hid_t file_space = H5Dget_space(did);
    hsize_t loc[1] = { dims[0]-1 };
    H5Sselect_elements( file_space, H5S_SELECT_SET, 1, &loc[0] );
    // Define the space in memory as a single int
    hsize_t onetime[1] = { 1 };
    hid_t memspace = H5Screate_simple(1, onetime, NULL);
    // Write the current timestep
    H5Dwrite( did, H5T_NATIVE_INT, memspace , file_space , transfer, &timestep );
    H5Sclose( memspace );
    H5Sclose(file_space);
    H5Dclose(did);
    
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    
    // Clear buffers
    data_uint  .resize(0);
    data_short .resize(0);
    data_double.resize(0);

}

