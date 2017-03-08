
#include <string>
#include <sstream>

#include "DiagnosticTrack.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, int speciesId ) :
IDs_done( params.restart ),
nDim_particle(params.nDim_particle)
{
    speciesId_ = speciesId;
    Species* species = patch->vecSpecies[speciesId_];
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE);
    
    ostringstream name("");
    name << "Tracking species '" << species->species_type << "'";
    
    // Get parameter "track_every" which describes a iteration selection
    timeSelection = new TimeSelection( PyTools::extract_py("track_every", "Species", speciesId_), name.str() );
    
    // Get parameter "track_flush_every" which decides the file flushing time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py("track_flush_every", "Species", speciesId_), name.str() );
    
    // Create the filename
    ostringstream hdf_filename("");
    hdf_filename << "TrackParticlesDisordered_" << species->species_type  << ".h5" ;
    filename = hdf_filename.str();
    
    // Create a list of the necessary datasets
    datasets.push_back( "Id" );
    datatypes.push_back( H5T_NATIVE_UINT64 );
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
    for (unsigned int idim=0 ; idim<nDim_particle ; idim++) {
        name.str("");
        name << "Position-" << idim;
        datasets.push_back( name.str() );
        datatypes.push_back( H5T_NATIVE_DOUBLE );
    }
    
}

DiagnosticTrack::~DiagnosticTrack()
{
    delete timeSelection;
    delete flush_timeSelection;
    H5Pclose(transfer);
}


void DiagnosticTrack::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    
    if ( newfile ) {
        // Create HDF5 file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Create the dataset for the time array
        hsize_t ntimes[1] = { 0 };
        hsize_t maxtimes[1] = { H5S_UNLIMITED };
        hid_t file_space = H5Screate_simple(1, ntimes, maxtimes);
        hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_layout(plist, H5D_CHUNKED);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        hsize_t chunks = 1;
        H5Pset_chunk(plist, 1, &chunks);
        hid_t did = H5Dcreate(fileId_, "Times", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dclose(did);
        H5Pclose(plist);
        H5Sclose(file_space);
        
    }
    else {
        // Open the file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid);
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


void DiagnosticTrack::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // Set the IDs of the particles
    if( ! IDs_done ) {
        
        latest_Id = smpi->getRank() * 4294967296; // 2^32
        
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            unsigned int s = vecPatches(ipatch)->vecSpecies[speciesId_]->particles->size();
            for (unsigned int iPart=0; iPart<s; iPart++) {
                latest_Id++;
                vecPatches(ipatch)->vecSpecies[speciesId_]->particles->id(iPart) = latest_Id;
            }
        }
        
        IDs_done = true;
    }
    
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    
}


bool DiagnosticTrack::prepare( int itime )
{
    return timeSelection->theTimeIsNow(itime);
}


void DiagnosticTrack::run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime )
{
    uint32_t nParticles_local = 0;
    uint64_t nParticles_global;
    unsigned int nPatches = vecPatches.size();
    
    hid_t current_group=0, plist=0, file_space=0, mem_space=0;
    #pragma omp master
    {
        // Obtain the particle partition of all the patches in this MPI
        patch_start.resize( vecPatches.size() );
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            patch_start[ipatch] = nParticles_local;
            nParticles_local += vecPatches(ipatch)->vecSpecies[speciesId_]->getNbrOfParticles();
        }
        
        // Specify the memory dataspace (the size of the local buffer)
        hsize_t count_[1] = {(hsize_t)nParticles_local};
        mem_space = H5Screate_simple(1, count_, NULL);
        
        // Get the number of particles for each MPI
        int sz = smpi->getSize(), rk = smpi->getRank();
        std::vector<uint32_t> all_nPart(sz);
        MPI_Allgather( &nParticles_local, 1, MPI_UNSIGNED, &all_nPart[0], 1, MPI_UNSIGNED, MPI_COMM_WORLD );
        
        // Calculate the cumulative sum, which is where this MPI will start writing
        uint64_t offset = 0;
        for (int irk=0; irk<rk; irk++) offset += all_nPart[irk];
        
        // Continue to cumulate in order to get the total number of particles
        nParticles_global = offset;
        for (int irk=rk; irk<sz; irk++) nParticles_global += all_nPart[irk];
        
        // Make a new group for this iteration
        ostringstream t("");
        t << setfill('0') << setw(10) << itime;
        current_group = H5::group( fileId_, t.str().c_str() );
        
        // Set the dataset parameters
        plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump
        H5Pset_layout(plist, H5D_CHUNKED);
        
        // Set the chunk size
        unsigned int maximum_chunk_size = 100000000;
        unsigned int number_of_chunks = nParticles_global/maximum_chunk_size;
        if( nParticles_global%maximum_chunk_size != 0 ) number_of_chunks++;
        unsigned int chunk_size = nParticles_global/number_of_chunks;
        if( nParticles_global%number_of_chunks != 0 ) chunk_size++;
        hsize_t chunk_dims = chunk_size;
        H5Pset_chunk(plist, 1, &chunk_dims);
        
        // Define maximum size
        hsize_t dims = nParticles_global;
        file_space = H5Screate_simple(1, &dims, NULL);   
        
        // Select locations that this proc will write
        if(nParticles_local>0) {
            hsize_t start=offset, count=1, block=nParticles_local;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, &start, NULL, &count, &block );
        } else {
            H5Sselect_none(file_space);
        }
    }
    
    #pragma omp master
    {
        // Create the "latest_IDs" dataset
        // Create file space and select one element for each proc
        hsize_t numel = smpi->getSize();
        hid_t filespace = H5Screate_simple(1, &numel, NULL);
        hsize_t offset=smpi->getRank(), count=1;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, NULL);
        // Create dataset
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hid_t dset_id  = H5Dcreate( current_group, "latest_IDs", H5T_NATIVE_UINT64, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
        // Create memory space
        hsize_t size_in_memory = 1;
        hid_t memspace  = H5Screate_simple(1, &size_in_memory, NULL);
        // Parallel write
        hid_t write_plist = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
        H5Dwrite( dset_id, H5T_NATIVE_UINT64, memspace, filespace, write_plist, &latest_Id );
        // Close all
        H5Pclose( write_plist );
        H5Pclose( plist_id );
        H5Dclose( dset_id );
        H5Sclose( filespace );
        H5Sclose( memspace );
    }
    
    // For each dataset
    for( unsigned int idset=0; idset<datasets.size(); idset++) {
        
        // Fill the buffer for the current patch
        unsigned int patch_nParticles, i, j;
        Particles* particles;
        if( idset == 0 ) { // Id
            #pragma omp master
            data_uint64.resize( nParticles_local, 0 );
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
                particles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                patch_nParticles = particles->size();
                i=0;
                j=patch_start[ipatch];
                while( i<patch_nParticles ) {
                    data_uint64[j] = particles->id(i);
                    i++; j++;
                }
            }
        } else if( idset == 1 ) { // Charge
            #pragma omp master
            data_short.resize( nParticles_local );
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
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
            #pragma omp master
            data_double.resize( nParticles_local );
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
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
            #pragma omp master
            data_double.resize( nParticles_local );
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
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
            #pragma omp master
            data_double.resize( nParticles_local );
            #pragma omp barrier
            #pragma omp for schedule(runtime)
            for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
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
        
        #pragma omp master
        {
            // Create the dataset
            hid_t did = H5Dcreate(current_group, datasets[idset].c_str(), datatypes[idset], file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
            
            // Write
            if( idset == 0 ) {
                H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_uint64[0] );
            } else if( idset == 1 ) {
                H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_short[0] );
            } else {
                H5Dwrite( did, datatypes[idset], mem_space , file_space , transfer, &data_double[0] );
            }
            
            H5Dclose(did);
        }
    }
    
    #pragma omp master
    {
        H5Pclose(plist);
        H5Sclose(file_space);
        H5Sclose( mem_space );
        H5Gclose( current_group );
        
        // Update the size of the times dataset for the new iteration
        hid_t did = H5Dopen( fileId_, "Times", H5P_DEFAULT );
        hid_t space = H5Dget_space( did );
        hssize_t ntimes = H5Sget_simple_extent_npoints( space );
        hsize_t tdims = ntimes+1;
        H5Dset_extent(did, &tdims);
        H5Sclose(space);
        // Select only the last element of the array
        space = H5Dget_space( did );
        hsize_t loc = ntimes;
        H5Sselect_elements( space, H5S_SELECT_SET, 1, &loc );
        // Define the space in memory as a single int
        hsize_t onetime = 1;
        hid_t memspace = H5Screate_simple(1, &onetime, NULL);
        // Write the current iteration
        H5Dwrite( did, H5T_NATIVE_INT, memspace , space , transfer, &itime );
        H5Sclose( memspace );
        H5Sclose(space);
        H5Dclose(did);
        
        if( flush_timeSelection->theTimeIsNow(itime) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        
        // Clear buffers
        data_uint64.resize(0);
        data_short .resize(0);
        data_double.resize(0);
    }
    #pragma omp barrier
}

