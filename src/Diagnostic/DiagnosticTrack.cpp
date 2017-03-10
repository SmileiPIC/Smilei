
#include <string>
#include <sstream>

#include "DiagnosticTrack.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, Patch* patch, unsigned int speciesId, OpenPMDparams& oPMD ) :
    Diagnostic(oPMD),
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
        
        // Attributes for openPMD
        openPMD->writeRootAttributes( fileId_, "no_meshes", "particles/" );
        
        // Create "data" group for openPMD compatibility
        data_group_id = H5::group(fileId_, "data");
        
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
        H5Gclose( data_group_id );
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
    uint64_t nParticles_global = 0;
    string xyz = "xyz";
    
    hid_t momentum_group=0, position_group=0, iteration_group=0, particles_group=0, species_group=0;
    hid_t plist=0, file_space=0, mem_space=0;
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
        iteration_group = H5::group( data_group_id, t.str().c_str() );
        particles_group = H5::group( iteration_group, "particles" );
        species_group = H5::group( particles_group, vecPatches(0)->vecSpecies[speciesId_]->species_type.c_str() );
        
        // Add openPMD attributes ( "basePath" )
        openPMD->writeBasePathAttributes( iteration_group, itime );
        // Add openPMD attributes ( "particles" )
        openPMD->writeParticlesAttributes( particles_group );
        // Add openPMD attributes ( path of a given species )
        openPMD->writeSpeciesAttributes( species_group );
        
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
        hid_t dset_id  = H5Dcreate( iteration_group, "latest_IDs", H5T_NATIVE_UINT64, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
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
    
    // Id
    #pragma omp master
    data_uint64.resize( nParticles_local, 1 );
    #pragma omp barrier
    fill_buffer(vecPatches, 0, data_uint64);
    #pragma omp master
    {
        write_scalar( species_group, "id", data_uint64[0], H5T_NATIVE_UINT64, file_space, mem_space, plist, 0 );
        data_uint64.resize(0);
    }
    
    // Charge
    #pragma omp master
    data_short.resize( nParticles_local, 0 );
    #pragma omp barrier
    fill_buffer(vecPatches, 0, data_short);
    #pragma omp master
    {
        write_scalar( species_group, "charge", data_short[0], H5T_NATIVE_SHORT, file_space, mem_space, plist, 7 );
        data_short.resize(0);
    }
    
    // Weight
    #pragma omp master
    data_double.resize( nParticles_local, 0 );
    #pragma omp barrier
    fill_buffer(vecPatches, nDim_particle+3, data_double);
    #pragma omp master
    write_scalar( species_group, "weight", data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, 4 );
    
    // Momentum
    #pragma omp master
    {
        momentum_group = H5::group(species_group, "momentum");
        openPMD->writeRecordAttributes( momentum_group, 6 );
    }
    for( unsigned int idim=0; idim<3; idim++ ) {
        #pragma omp barrier
        fill_buffer(vecPatches, nDim_particle+idim, data_double);
        #pragma omp master
        write_component( momentum_group, xyz.substr(idim,1).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist );
    }
    #pragma omp master
    H5Gclose( momentum_group );
    
    // Position
    #pragma omp master
    {
        position_group = H5::group(species_group, "position");
        openPMD->writeRecordAttributes( position_group, 5 );
    }
    for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
        #pragma omp barrier
        fill_buffer(vecPatches, idim, data_double);
        #pragma omp master
        write_component( position_group, xyz.substr(idim,1).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist );
    }
    #pragma omp master
    H5Gclose( position_group );
    
    #pragma omp master
    {
        
    // PositionOffset (for OpenPMD)
        hid_t positionoffset_group = H5::group(species_group, "positionOffset");
        openPMD->writeRecordAttributes( positionoffset_group, 5 );
        vector<uint64_t> np = {nParticles_global};
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            hid_t xyz_group = H5::group(positionoffset_group, xyz.substr(idim,1));
            openPMD->writeComponentAttributes( xyz_group );
            H5::attr( xyz_group, "value", 0. );
            H5::attr( xyz_group, "shape", np, H5T_NATIVE_UINT64 );
            H5Gclose( xyz_group );
        }
        H5Gclose( positionoffset_group );
        
    // Close and flush
        data_double.resize(0);
        
        H5Pclose(plist);
        H5Sclose(file_space);
        H5Sclose( mem_space );
        H5Gclose( species_group );
        H5Gclose( particles_group );
        H5Gclose( iteration_group );
        
        if( flush_timeSelection->theTimeIsNow(itime) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
    #pragma omp barrier
}


template<typename T>
void DiagnosticTrack::fill_buffer(VectorPatch& vecPatches, unsigned int iprop, vector<T>& buffer)
{
    unsigned int patch_nParticles, i, j, nPatches=vecPatches.size();
    vector<T>* property = NULL;
    #pragma omp for schedule(runtime)
    for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
        patch_nParticles = vecPatches(ipatch)->vecSpecies[speciesId_]->particles->size();
        vecPatches(ipatch)->vecSpecies[speciesId_]->particles->getProperty(iprop, property);
        i=0;
        j=patch_start[ipatch];
        while( i<patch_nParticles ) {
            buffer[j] = (*property)[i];
            i++; j++;
        }
    }
}


template<typename T>
void DiagnosticTrack::write_scalar( hid_t location, string name, T& buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist, unsigned int unit_type )
{
    hid_t did = H5Dcreate(location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dwrite( did, dtype, mem_space , file_space , transfer, &buffer );
    openPMD->writeRecordAttributes( did, unit_type );
    openPMD->writeComponentAttributes( did );
    H5Dclose(did);
}

template<typename T>
void DiagnosticTrack::write_component( hid_t location, string name, T& buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist )
{
    hid_t did = H5Dcreate(location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    H5Dwrite( did, dtype, mem_space , file_space , transfer, &buffer );
    openPMD->writeComponentAttributes( did );
    H5Dclose(did);
}
