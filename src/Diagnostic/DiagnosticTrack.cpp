#include "PyTools.h"

#include <string>
#include <sstream>

#include "ParticleData.h"
#include "DiagnosticTrack.h"
#include "VectorPatch.h"
#include "Params.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, unsigned int iDiagTrackParticles, unsigned int idiag, OpenPMDparams& oPMD ) :
    Diagnostic(oPMD),
    IDs_done( params.restart ),
    nDim_particle(params.nDim_particle)
{
    // Extract the species
    string species_name;
    if( !PyTools::extract("species",species_name,"DiagTrackParticles",iDiagTrackParticles) )
        ERROR("DiagTrackParticles #" << iDiagTrackParticles << " requires an argument `species`");
    vector<string> species_names = {species_name};
    vector<unsigned int> species_ids = Params::FindSpecies(vecPatches(0)->vecSpecies, species_names);
    if( species_ids.size() > 1 )
        ERROR("DiagTrackParticles #" << iDiagTrackParticles << " corresponds to more than 1 species");
    if( species_ids.size() < 1 )
        ERROR("DiagTrackParticles #" << iDiagTrackParticles << " does not correspond to any existing species");
    speciesId_ = species_ids[0];
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE);

    ostringstream name("");
    name << "Tracking species '" << species_name << "'";
    
    // Get parameter "every" which describes an iteration selection
    timeSelection = new TimeSelection( PyTools::extract_py("every", "DiagTrackParticles", iDiagTrackParticles), name.str() );
    
    // Get parameter "flush_every" which decides the file flushing time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py("flush_every", "DiagTrackParticles", iDiagTrackParticles), name.str() );
    
    // Inform each patch about this diag
    for( unsigned int ipatch=0; ipatch<vecPatches.size(); ipatch++ ) {
        vecPatches(ipatch)->vecSpecies[speciesId_]->tracking_diagnostic = idiag;
    }
    
    // Get parameter "filter" which gives a python function to select particles
    filter = PyTools::extract_py("filter", "DiagTrackParticles", iDiagTrackParticles);
    has_filter = (filter != Py_None);
    if( has_filter ) {
#ifdef SMILEI_USE_NUMPY
        PyTools::setIteration( 0 );
        // Test the filter with temporary, "fake" particles
        name << " filter:";
        bool * dummy = NULL;
        ParticleData test( nDim_particle, filter, name.str(), dummy );
#else
        ERROR(name.str() << " with a filter requires the numpy package");
#endif
    }

    // Create the filename
    ostringstream hdf_filename("");
    hdf_filename << "TrackParticlesDisordered_" << species_name  << ".h5" ;
    filename = hdf_filename.str();

}

DiagnosticTrack::~DiagnosticTrack()
{
    delete timeSelection;
    delete flush_timeSelection;
    H5Pclose(transfer);
    Py_DECREF(filter);
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

        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++)
            setIDs(vecPatches(ipatch));

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


void DiagnosticTrack::run( SmileiMPI* smpi, VectorPatch& vecPatches, int itime, SimWindow* simWindow )
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
        
        if( has_filter ) {
        
#ifdef SMILEI_USE_NUMPY
            // Set a python variable "Main.iteration" to itime so that it can be accessed in the filter
            PyTools::setIteration( itime );
            
            patch_selection.resize( vecPatches.size() );
            PyArrayObject *ret;
            ParticleData particleData(0);
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                Particles * p = vecPatches(ipatch)->vecSpecies[speciesId_]->particles;
                unsigned int npart = p->size();
                // Expose particle data as numpy arrays
                particleData.resize( npart );
                particleData.set( p );
                // run the filter function
                ret = (PyArrayObject*)PyObject_CallFunctionObjArgs(filter, particleData.get(), NULL);
                particleData.clear();
                
                // Loop the return value and store the particle IDs
                bool* arr = (bool*) PyArray_GETPTR1( ret, 0 );
                patch_selection[ipatch].resize(0);
                for(unsigned int i=0; i<npart; i++) {
                    if( arr[i] ) {
                        patch_selection[ipatch].push_back( i );
                        // If particle not tracked before (ID==0), then set its ID
                        if( p->id(i) == 0 ) p->id(i) = ++latest_Id;
                    }
                }
                patch_start[ipatch] = nParticles_local;
                nParticles_local += patch_selection[ipatch].size();
                
                Py_DECREF(ret);
            }
#endif
        
        } else {
            for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
                patch_start[ipatch] = nParticles_local;
                nParticles_local += vecPatches(ipatch)->vecSpecies[speciesId_]->getNbrOfParticles();
            }
        }
        
        // Specify the memory dataspace (the size of the local buffer)
        hsize_t count_[1] = {(hsize_t)nParticles_local};
        mem_space = H5Screate_simple(1, count_, NULL);
        
        // Get the number of offset for this MPI rank
        uint64_t np_local = nParticles_local, offset;
        MPI_Scan( &np_local, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
        nParticles_global = offset;
        offset -= np_local;
        MPI_Bcast( &nParticles_global, 1, MPI_UNSIGNED_LONG_LONG, smpi->getSize()-1, MPI_COMM_WORLD );
        
        // Make a new group for this iteration
        ostringstream t("");
        t << setfill('0') << setw(10) << itime;
        iteration_group = H5::group( data_group_id, t.str().c_str() );
        particles_group = H5::group( iteration_group, "particles" );
        species_group = H5::group( particles_group, vecPatches(0)->vecSpecies[speciesId_]->name.c_str() );
        
        // Add openPMD attributes ( "basePath" )
        openPMD->writeBasePathAttributes( iteration_group, itime );
        // Add openPMD attributes ( "particles" )
        openPMD->writeParticlesAttributes( particles_group );
        // Add openPMD attributes ( path of a given species )
        openPMD->writeSpeciesAttributes( species_group );

        // Set the dataset parameters
        plist = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_alloc_time(plist, H5D_ALLOC_TIME_EARLY); // necessary for collective dump

        if( nParticles_global>0 ){
            H5Pset_layout(plist, H5D_CHUNKED);

            // Set the chunk size
            unsigned int maximum_chunk_size = 100000000;
            unsigned int number_of_chunks = nParticles_global/maximum_chunk_size;
            if( nParticles_global%maximum_chunk_size != 0 ) number_of_chunks++;
            if( number_of_chunks==0 ) number_of_chunks = 1;
            unsigned int chunk_size = nParticles_global/number_of_chunks;
            if( nParticles_global%number_of_chunks != 0 ) chunk_size++;
            hsize_t chunk_dims = chunk_size;
            H5Pset_chunk(plist, 1, &chunk_dims);
        }

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

        // Create the "latest_IDs" dataset
        // Create file space and select one element for each proc
        hsize_t numel = smpi->getSize();
        hid_t filespace = H5Screate_simple(1, &numel, NULL);
        hsize_t offset_ = smpi->getRank(), count=1;
        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset_, NULL, &count, NULL);
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
        write_scalar( species_group, "id", data_uint64[0], H5T_NATIVE_UINT64, file_space, mem_space, plist, SMILEI_UNIT_NONE, nParticles_global );
        data_uint64.resize(0);
    }

    // Charge
    #pragma omp master
    data_short.resize( nParticles_local, 0 );
    #pragma omp barrier
    fill_buffer(vecPatches, 0, data_short);
    #pragma omp master
    {
        write_scalar( species_group, "charge", data_short[0], H5T_NATIVE_SHORT, file_space, mem_space, plist, SMILEI_UNIT_CHARGE, nParticles_global );
        data_short.resize(0);
    }

    // Weight
    #pragma omp master
    data_double.resize( nParticles_local, 0 );
    #pragma omp barrier
    fill_buffer(vecPatches, nDim_particle+3, data_double);
    #pragma omp master
    write_scalar( species_group, "weight", data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_DENSITY, nParticles_global );

    // Momentum
    #pragma omp master
    {
        momentum_group = H5::group(species_group, "momentum");
        openPMD->writeRecordAttributes( momentum_group, SMILEI_UNIT_MOMENTUM );
    }
    for( unsigned int idim=0; idim<3; idim++ ) {
        #pragma omp barrier
        fill_buffer(vecPatches, nDim_particle+idim, data_double);
        #pragma omp master
        write_component( momentum_group, xyz.substr(idim,1).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_MOMENTUM, nParticles_global );
    }
    #pragma omp master
    H5Gclose( momentum_group );

    // Position
    #pragma omp master
    {
        position_group = H5::group(species_group, "position");
        openPMD->writeRecordAttributes( position_group, SMILEI_UNIT_POSITION );
    }
    for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
        #pragma omp barrier
        fill_buffer(vecPatches, idim, data_double);
        #pragma omp master
        write_component( position_group, xyz.substr(idim,1).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_POSITION, nParticles_global );
    }
    #pragma omp master
    H5Gclose( position_group );

    // If the quantum parameter is available
    if (vecPatches(0)->vecSpecies[speciesId_]->particles->isQuantumParameter)
    {
        // Chi - quantum parameter
        #pragma omp master
        data_double.resize( nParticles_local, 0 );
        #pragma omp barrier
// Position old exists in this case
#ifdef  __DEBUG
        fill_buffer(vecPatches, nDim_particle+3+3+1, data_double);
// Else, position old exists in this case
#else
        fill_buffer(vecPatches, nDim_particle+3+1, data_double);
#endif
        #pragma omp master
        {
            write_scalar( species_group, "chi", data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_NONE, nParticles_global );
        }
    }

    #pragma omp master
    {

    // PositionOffset (for OpenPMD)
        hid_t positionoffset_group = H5::group(species_group, "positionOffset");
        openPMD->writeRecordAttributes( positionoffset_group, SMILEI_UNIT_POSITION );
        vector<uint64_t> np = {nParticles_global};
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            hid_t xyz_group = H5::group(positionoffset_group, xyz.substr(idim,1));
            openPMD->writeComponentAttributes( xyz_group, SMILEI_UNIT_POSITION );
            H5::attr( xyz_group, "value", 0. );
            H5::attr( xyz_group, "shape", np, H5T_NATIVE_UINT64 );
            H5Gclose( xyz_group );
        }
        H5Gclose( positionoffset_group );

    // Close and flush
        data_double.resize(0);
        patch_selection.resize(0);

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


void DiagnosticTrack::setIDs(Patch * patch)
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) return;
    unsigned int s = patch->vecSpecies[speciesId_]->particles->size();
    for (unsigned int iPart=0; iPart<s; iPart++)
        patch->vecSpecies[speciesId_]->particles->id(iPart) = ++latest_Id;
}


void DiagnosticTrack::setIDs(Particles& particles)
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) return;
    unsigned int s = particles.size(), id;
    #pragma omp critical
    {
        for (unsigned int iPart=0; iPart<s; iPart++) {
            id = ++latest_Id;
            particles.id(iPart) = id;
        }
    }
}


template<typename T>
void DiagnosticTrack::fill_buffer(VectorPatch& vecPatches, unsigned int iprop, vector<T>& buffer)
{
    unsigned int patch_nParticles, i, j, nPatches=vecPatches.size();
    vector<T>* property = NULL;

    if( has_filter ) {
        #pragma omp for schedule(runtime)
        for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++) {
            patch_nParticles = patch_selection[ipatch].size();
            vecPatches(ipatch)->vecSpecies[speciesId_]->particles->getProperty(iprop, property);
            i=0;
            j=patch_start[ipatch];
            while( i<patch_nParticles ) {
                buffer[j] = (*property)[patch_selection[ipatch][i]];
                i++; j++;
            }
        }
    } else {
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
}


template<typename T>
void DiagnosticTrack::write_scalar( hid_t location, string name, T& buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist, unsigned int unit_type, unsigned int npart_global )
{
    hid_t did = H5Dcreate(location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if( npart_global>0 ) H5Dwrite( did, dtype, mem_space , file_space , transfer, &buffer );
    openPMD->writeRecordAttributes( did, unit_type );
    openPMD->writeComponentAttributes( did, unit_type );
    H5Dclose(did);
}

template<typename T>
void DiagnosticTrack::write_component( hid_t location, string name, T& buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist, unsigned int unit_type, unsigned int npart_global )
{
    hid_t did = H5Dcreate(location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
    if( npart_global>0 ) H5Dwrite( did, dtype, mem_space , file_space , transfer, &buffer );
    openPMD->writeComponentAttributes( did, unit_type );
    H5Dclose(did);
}
