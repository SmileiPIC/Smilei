#include "PyTools.h"

#include <string>
#include <sstream>

#include "ParticleData.h"
#include "PeekAtSpecies.h"
#include "DiagnosticTrack.h"
#include "VectorPatch.h"
#include "Params.h"

using namespace std;

DiagnosticTrack::DiagnosticTrack( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int iDiagTrackParticles, unsigned int idiag, OpenPMDparams &oPMD ) :
    Diagnostic( &oPMD, "DiagTrackParticles", iDiagTrackParticles ),
    IDs_done( params.restart ),
    nDim_particle( params.nDim_particle )
{

    // Extract the species
    string species_name;
    PyTools::extract( "species", species_name, "DiagTrackParticles", iDiagTrackParticles );
    vector<string> species_names = {species_name};
    vector<unsigned int> species_ids = Params::FindSpecies( vecPatches( 0 )->vecSpecies, species_names );
    if( species_ids.size() > 1 ) {
        ERROR( "DiagTrackParticles #" << iDiagTrackParticles << " corresponds to more than 1 species" );
    }
    if( species_ids.size() < 1 ) {
        ERROR( "DiagTrackParticles #" << iDiagTrackParticles << " does not correspond to any existing species" );
    }
    speciesId_ = species_ids[0];
    
    // Define the transfer type (collective is faster than independent)
    transfer = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( transfer, H5FD_MPIO_COLLECTIVE );
    
    ostringstream name( "" );
    name << "Tracking species '" << species_name << "'";
    
    // Get parameter "every" which describes an iteration selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagTrackParticles", iDiagTrackParticles ), name.str() );
    
    // Get parameter "flush_every" which decides the file flushing time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py( "flush_every", "DiagTrackParticles", iDiagTrackParticles ), name.str() );
    
    // Inform each patch about this diag
    for( unsigned int ipatch=0; ipatch<vecPatches.size(); ipatch++ ) {
        vecPatches( ipatch )->vecSpecies[speciesId_]->tracking_diagnostic = idiag;
    }
    
    // Get parameter "filter" which gives a python function to select particles
    filter = PyTools::extract_py( "filter", "DiagTrackParticles", iDiagTrackParticles );
    has_filter = ( filter != Py_None );
    if( has_filter ) {
#ifdef SMILEI_USE_NUMPY
        PyTools::setIteration( 0 );
        // Test the filter with temporary, "fake" particles
        name << " filter:";
        bool *dummy = NULL;
        ParticleData test( nDim_particle, filter, name.str(), dummy );
#else
        ERROR( name.str() << " with a filter requires the numpy package" );
#endif
    }
    
    // Get the parameter "attributes": a list of attribute name that must be written
    vector<string> attributes( 0 );
    if( !PyTools::extractV( "attributes", attributes, "DiagTrackParticles", iDiagTrackParticles ) ) {
        ERROR( "DiagTrackParticles #" << iDiagTrackParticles << ": argument `attribute` must be a list of strings" );
    }
    if( attributes.size() == 0 ) {
        ERROR( "DiagTrackParticles #" << iDiagTrackParticles << ": argument `attribute` must have at least one element" );
    }
    ostringstream attr_list( "" );
    attr_list << "id";
    write_position.resize( 3, false );
    write_momentum.resize( 3, false );
    write_charge = false;
    write_weight = false;
    write_chi    = false;
    write_E.resize( 3, false );
    write_B.resize( 3, false );
    interpolate = false;
    for( unsigned int i=0; i<attributes.size(); i++ ) {
        if( attributes[i] == "x" ) {
            write_position[0] = true;
        } else if( attributes[i] == "y" ) {
            if( nDim_particle>1 ) {
                write_position[1] = true;
            } else {
                continue;
            }
        } else if( attributes[i] == "z" ) {
            if( nDim_particle>2 ) {
                write_position[2] = true;
            } else {
                continue;
            }
        } else if( attributes[i] == "px" ) {
            write_momentum[0] = true;
        } else if( attributes[i] == "py" ) {
            write_momentum[1] = true;
        } else if( attributes[i] == "pz" ) {
            write_momentum[2] = true;
        } else if( attributes[i] == "charge" || attributes[i] == "q" ) {
            write_charge      = true;
        } else if( attributes[i] == "weight" || attributes[i] == "w" ) {
            write_weight      = true;
        } else if( attributes[i] == "chi" ) {
            write_chi         = true;
        } else if( attributes[i] == "Ex" ) {
            write_E[0]        = true;
            interpolate = true;
        } else if( attributes[i] == "Ey" ) {
            write_E[1]        = true;
            interpolate = true;
        } else if( attributes[i] == "Ez" ) {
            write_E[2]        = true;
            interpolate = true;
        } else if( attributes[i] == "Bx" ) {
            write_B[0]        = true;
            interpolate = true;
        } else if( attributes[i] == "By" ) {
            write_B[1]        = true;
            interpolate = true;
        } else if( attributes[i] == "Bz" ) {
            write_B[2]        = true;
            interpolate = true;
        } else {
            ERROR( "DiagTrackParticles #" << iDiagTrackParticles << ": attribute `" << attributes[i] << "` unknown" );
        }
        attr_list << "," << attributes[i];
    }
    write_any_position = write_position[0] || write_position[1] || write_position[2];
    write_any_momentum = write_momentum[0] || write_momentum[1] || write_momentum[2];
    write_any_E = write_E[0] || write_E[1] || write_E[2];
    write_any_B = write_B[0] || write_B[1] || write_B[2];
    if( write_chi && ! vecPatches( 0 )->vecSpecies[speciesId_]->particles->isQuantumParameter ) {
        ERROR( "DiagTrackParticles #" << iDiagTrackParticles << ": attribute `chi` not available for this species" );
    }
    
    // Create the filename
    ostringstream hdf_filename( "" );
    hdf_filename << "TrackParticlesDisordered_" << species_name  << ".h5" ;
    filename = hdf_filename.str();
    
    // Print some info
    if( smpi->isMaster() ) {
        MESSAGE( 1, "Created TrackParticles #" << iDiagTrackParticles << ": species " << species_name );
        MESSAGE( 2, attr_list.str() );
    }
    
    // Obtain the approximate number of particles in the species
    if( params.print_expected_disk_usage ) {
        PeekAtSpecies peek( params, speciesId_ );
        npart_total = peek.totalNumberofParticles();
    } else {
        npart_total = 0;
    }
}

DiagnosticTrack::~DiagnosticTrack()
{
    delete timeSelection;
    delete flush_timeSelection;
    H5Pclose( transfer );
    Py_DECREF( filter );
}


void DiagnosticTrack::openFile( Params &params, SmileiMPI *smpi, bool newfile )
{

    if( newfile ) {
        // Create HDF5 file
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_ = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid );
        H5Pclose( pid );
        
        H5::attr( fileId_, "name", diag_name_ );
        
        // Attributes for openPMD
        openPMD_->writeRootAttributes( fileId_, "no_meshes", "particles/" );
        
        // Create "data" group for openPMD compatibility
        data_group_id = H5::group( fileId_, "data" );
        
    } else {
        // Open the file
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose( pid );
    }
    
}


void DiagnosticTrack::closeFile()
{
    if( fileId_>0 ) {
        H5Gclose( data_group_id );
        H5Fclose( fileId_ );
        fileId_=0;
    }
}


void DiagnosticTrack::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // Set the IDs of the particles
    if( ! IDs_done ) {
        latest_Id = smpi->getRank() * 4294967296; // 2^32
        
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
            setIDs( vecPatches( ipatch ) );
        }
        
        IDs_done = true;
    }
    
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    
}


bool DiagnosticTrack::prepare( int itime )
{
    return timeSelection->theTimeIsNow( itime );
}


void DiagnosticTrack::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers )
{
    uint64_t nParticles_global = 0;
    string xyz = "xyz";
    
    hid_t momentum_group=0, position_group=0, iteration_group=0, particles_group=0, species_group=0;
    hid_t plist=0, file_space=0, mem_space=0;
    #pragma omp master
    {
        // Obtain the particle partition of all the patches in this MPI
        nParticles_local = 0;
        patch_start.resize( vecPatches.size() );
        
        if( has_filter ) {
        
#ifdef SMILEI_USE_NUMPY
            // Set a python variable "Main.iteration" to itime so that it can be accessed in the filter
            PyTools::setIteration( itime );
            
            patch_selection.resize( vecPatches.size() );
            PyArrayObject *ret;
            ParticleData particleData( 0 );
            for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
                patch_selection[ipatch].resize( 0 );
                Particles *p = vecPatches( ipatch )->vecSpecies[speciesId_]->particles;
                unsigned int npart = p->size();
                if( npart > 0 ) {
                    // Expose particle data as numpy arrays
                    particleData.resize( npart );
                    particleData.set( p );
                    // run the filter function
                    ret = ( PyArrayObject * )PyObject_CallFunctionObjArgs( filter, particleData.get(), NULL );
                    PyTools::checkPyError();
                    particleData.clear();
                    if( ret == NULL ) {
                        ERROR( "A DiagTrackParticles filter has not provided a correct result" );
                    }
                    // Loop the return value and store the particle IDs
                    bool *arr = ( bool * ) PyArray_GETPTR1( ret, 0 );
                    for( unsigned int i=0; i<npart; i++ ) {
                        if( arr[i] ) {
                            patch_selection[ipatch].push_back( i );
                            // If particle not tracked before (ID==0), then set its ID
                            if( p->id( i ) == 0 ) {
                                p->id( i ) = ++latest_Id;
                            }
                        }
                    }
                    Py_DECREF( ret );
                }
                patch_start[ipatch] = nParticles_local;
                nParticles_local += patch_selection[ipatch].size();
            }
#endif
            
        } else {
            for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
                patch_start[ipatch] = nParticles_local;
                nParticles_local += vecPatches( ipatch )->vecSpecies[speciesId_]->getNbrOfParticles();
            }
        }
        
        // Specify the memory dataspace (the size of the local buffer)
        hsize_t count_ = nParticles_local;
        mem_space = H5Screate_simple( 1, &count_, NULL );
        
        // Get the number of offset for this MPI rank
        uint64_t np_local = nParticles_local, offset;
        MPI_Scan( &np_local, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD );
        nParticles_global = offset;
        offset -= np_local;
        MPI_Bcast( &nParticles_global, 1, MPI_UNSIGNED_LONG_LONG, smpi->getSize()-1, MPI_COMM_WORLD );
        
        // Make a new group for this iteration
        ostringstream t( "" );
        t << setfill( '0' ) << setw( 10 ) << itime;
        iteration_group = H5::group( data_group_id, t.str().c_str() );
        particles_group = H5::group( iteration_group, "particles" );
        species_group = H5::group( particles_group, vecPatches( 0 )->vecSpecies[speciesId_]->name_.c_str() );
        
        // Add openPMD attributes ( "basePath" )
        openPMD_->writeBasePathAttributes( iteration_group, itime );
        // Add openPMD attributes ( "particles" )
        openPMD_->writeParticlesAttributes( particles_group );
        // Add openPMD attributes ( path of a given species )
        openPMD_->writeSpeciesAttributes( species_group );
        
        // Write x_moved
        H5::attr( iteration_group, "x_moved", simWindow ? simWindow->getXmoved() : 0. );
        
        // Set the dataset parameters
        plist = H5Pcreate( H5P_DATASET_CREATE );
        H5Pset_alloc_time( plist, H5D_ALLOC_TIME_EARLY ); // necessary for collective dump
        
        if( nParticles_global>0 ) {
            // Set the chunk size
            unsigned int maximum_chunk_size = 100000000;
            unsigned int number_of_chunks = nParticles_global/maximum_chunk_size;
            if( nParticles_global%maximum_chunk_size != 0 ) {
                number_of_chunks++;
            }
            if( number_of_chunks==0 ) {
                number_of_chunks = 1;
            }
            unsigned int chunk_size = nParticles_global/number_of_chunks;
            if( nParticles_global%number_of_chunks != 0 ) {
                chunk_size++;
            }
            hsize_t chunk_dims = chunk_size;
            if( number_of_chunks > 1 ) {
                H5Pset_layout( plist, H5D_CHUNKED );
                H5Pset_chunk( plist, 1, &chunk_dims );
            }
        }
        
        // Define maximum size
        hsize_t dims = nParticles_global;
        file_space = H5Screate_simple( 1, &dims, NULL );
        
        // Select locations that this proc will write
        if( nParticles_local>0 ) {
            hsize_t start=offset, count=1, block=nParticles_local;
            H5Sselect_hyperslab( file_space, H5S_SELECT_SET, &start, NULL, &count, &block );
        } else {
            H5Sselect_none( file_space );
        }
        
        // Create the "latest_IDs" dataset
        // Create file space and select one element for each proc
        hsize_t numel = smpi->getSize();
        hid_t filespace = H5Screate_simple( 1, &numel, NULL );
        hsize_t offset_ = smpi->getRank(), count=1;
        H5Sselect_hyperslab( filespace, H5S_SELECT_SET, &offset_, NULL, &count, NULL );
        // Create dataset
        hid_t plist_id = H5Pcreate( H5P_DATASET_CREATE );
        hid_t dset_id  = H5Dcreate( iteration_group, "latest_IDs", H5T_NATIVE_UINT64, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
        // Create memory space
        hsize_t size_in_memory = 1;
        hid_t memspace  = H5Screate_simple( 1, &size_in_memory, NULL );
        // Parallel write
        hid_t write_plist = H5Pcreate( H5P_DATASET_XFER );
        H5Pset_dxpl_mpio( write_plist, H5FD_MPIO_COLLECTIVE );
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
    fill_buffer( vecPatches, 0, data_uint64 );
    #pragma omp master
    {
        write_scalar( species_group, "id", data_uint64[0], H5T_NATIVE_UINT64, file_space, mem_space, plist, SMILEI_UNIT_NONE, nParticles_global );
        data_uint64.resize( 0 );
    }
    
    // Charge
    if( write_charge ) {
        #pragma omp master
        data_short.resize( nParticles_local, 0 );
        #pragma omp barrier
        fill_buffer( vecPatches, 0, data_short );
        #pragma omp master
        {
            write_scalar( species_group, "charge", data_short[0], H5T_NATIVE_SHORT, file_space, mem_space, plist, SMILEI_UNIT_CHARGE, nParticles_global );
            data_short.resize( 0 );
        }
    }
    
    #pragma omp master
    data_double.resize( nParticles_local, 0 );
    
    // Weight
    if( write_weight ) {
        #pragma omp barrier
        fill_buffer( vecPatches, nDim_particle+3, data_double );
        #pragma omp master
        write_scalar( species_group, "weight", data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_DENSITY, nParticles_global );
    }
    
    // Momentum
    if( write_any_momentum ) {
        #pragma omp master
        {
            momentum_group = H5::group( species_group, "momentum" );
            openPMD_->writeRecordAttributes( momentum_group, SMILEI_UNIT_MOMENTUM );
        }
        for( unsigned int idim=0; idim<3; idim++ ) {
            if( write_momentum[idim] ) {
                #pragma omp barrier
                fill_buffer( vecPatches, nDim_particle+idim, data_double );
                #pragma omp master
                {
                    // Multiply by the mass to obtain an actual momentum (except for photons (mass = 0))
                    if( vecPatches( 0 )->vecSpecies[speciesId_]->mass_ != 1. &&
                        vecPatches( 0 )->vecSpecies[speciesId_]->mass_ > 0) {
                        for( unsigned int ip=0; ip<nParticles_local; ip++ ) {
                            data_double[ip] *= vecPatches( 0 )->vecSpecies[speciesId_]->mass_;
                        }
                    }
                    write_component( momentum_group, xyz.substr( idim, 1 ).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_MOMENTUM, nParticles_global );
                }
            }
        }
        #pragma omp master
        H5Gclose( momentum_group );
    }
    
    // Position
    if( write_any_position ) {
        #pragma omp master
        {
            position_group = H5::group( species_group, "position" );
            openPMD_->writeRecordAttributes( position_group, SMILEI_UNIT_POSITION );
        }
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            if( write_position[idim] ) {
                #pragma omp barrier
                fill_buffer( vecPatches, idim, data_double );
                #pragma omp master
                write_component( position_group, xyz.substr( idim, 1 ).c_str(), data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_POSITION, nParticles_global );
            }
        }
        #pragma omp master
        H5Gclose( position_group );
    }
    
    // Chi - quantum parameter
    if( write_chi ) {
        #pragma omp barrier
// Position old exists in this case
#ifdef  __DEBUG
        fill_buffer( vecPatches, nDim_particle+3+3+1, data_double );
// Else, position old does not exist
#else
        fill_buffer( vecPatches, nDim_particle+3+1, data_double );
#endif
        #pragma omp master
        write_scalar( species_group, "chi", data_double[0], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_NONE, nParticles_global );
    }
    
    #pragma omp barrier
    
    // If field interpolation necessary
    if( interpolate ) {
    
    
        #pragma omp master
        data_double.resize( nParticles_local*6 );
        
        // Do the interpolation
        unsigned int nPatches=vecPatches.size();
        #pragma omp barrier
        
        if( has_filter ) {
            #pragma omp for schedule(static)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                vecPatches.species( ipatch, speciesId_ )->Interp->fieldsSelection(
                    vecPatches.emfields( ipatch ),
                    *( vecPatches.species( ipatch, speciesId_ )->particles ),
                    &data_double[patch_start[ipatch]],
                    ( int ) nParticles_local,
                    &patch_selection[ipatch]
                );
            }
        } else {
            #pragma omp for schedule(static)
            for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
                vecPatches.species( ipatch, speciesId_ )->Interp->fieldsSelection(
                    vecPatches.emfields( ipatch ),
                    *( vecPatches.species( ipatch, speciesId_ )->particles ),
                    &data_double[patch_start[ipatch]],
                    ( int ) nParticles_local,
                    NULL
                );
            }
        }
        #pragma omp barrier
        
        // Write out the fields
        #pragma omp master
        {
            if( write_any_E ) {
                hid_t Efield_group = H5::group( species_group, "E" );
                openPMD_->writeRecordAttributes( Efield_group, SMILEI_UNIT_EFIELD );
                for( unsigned int idim=0; idim<3; idim++ ) {
                    if( write_E[idim] ) {
                        write_component( Efield_group, xyz.substr( idim, 1 ).c_str(), data_double[idim*nParticles_local], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_EFIELD, nParticles_global );
                    }
                }
                H5Gclose( Efield_group );
            }
            
            if( write_any_B ) {
                hid_t Bfield_group = H5::group( species_group, "B" );
                openPMD_->writeRecordAttributes( Bfield_group, SMILEI_UNIT_BFIELD );
                for( unsigned int idim=0; idim<3; idim++ ) {
                    if( write_B[idim] ) {
                        write_component( Bfield_group, xyz.substr( idim, 1 ).c_str(), data_double[( 3+idim )*nParticles_local], H5T_NATIVE_DOUBLE, file_space, mem_space, plist, SMILEI_UNIT_BFIELD, nParticles_global );
                    }
                }
                H5Gclose( Bfield_group );
            }
        }
    } // END if interpolate
    
    #pragma omp master
    {
        data_double.resize( 0 );
        
        // PositionOffset (for OpenPMD)
        hid_t positionoffset_group = H5::group( species_group, "positionOffset" );
        openPMD_->writeRecordAttributes( positionoffset_group, SMILEI_UNIT_POSITION );
        vector<uint64_t> np = {nParticles_global};
        for( unsigned int idim=0; idim<nDim_particle; idim++ ) {
            hid_t xyz_group = H5::group( positionoffset_group, xyz.substr( idim, 1 ) );
            openPMD_->writeComponentAttributes( xyz_group, SMILEI_UNIT_POSITION );
            H5::attr( xyz_group, "value", 0. );
            H5::attr( xyz_group, "shape", np, H5T_NATIVE_UINT64 );
            H5Gclose( xyz_group );
        }
        H5Gclose( positionoffset_group );
        
        // Close and flush
        patch_selection.resize( 0 );
        
        H5Pclose( plist );
        H5Sclose( file_space );
        H5Sclose( mem_space );
        H5Gclose( species_group );
        H5Gclose( particles_group );
        H5Gclose( iteration_group );
        
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        }
    }
    #pragma omp barrier
}


void DiagnosticTrack::setIDs( Patch *patch )
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) {
        return;
    }
    unsigned int s = patch->vecSpecies[speciesId_]->particles->size();
    for( unsigned int iPart=0; iPart<s; iPart++ ) {
        patch->vecSpecies[speciesId_]->particles->id( iPart ) = ++latest_Id;
    }
}


void DiagnosticTrack::setIDs( Particles &particles )
{
    // If filter, IDs are set on-the-fly
    if( has_filter ) {
        return;
    }
    unsigned int s = particles.size();
    #pragma omp critical
    {
        for( unsigned int iPart=0; iPart<s; iPart++ ) {
            particles.id( iPart ) = ++latest_Id;
        }
    }
}


template<typename T>
void DiagnosticTrack::fill_buffer( VectorPatch &vecPatches, unsigned int iprop, vector<T> &buffer )
{
    unsigned int patch_nParticles, i, j, nPatches=vecPatches.size();
    vector<T> *property = NULL;
    
    if( has_filter ) {
        #pragma omp for schedule(runtime)
        for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
            patch_nParticles = patch_selection[ipatch].size();
            vecPatches( ipatch )->vecSpecies[speciesId_]->particles->getProperty( iprop, property );
            i=0;
            j=patch_start[ipatch];
            while( i<patch_nParticles ) {
                buffer[j] = ( *property )[patch_selection[ipatch][i]];
                i++;
                j++;
            }
        }
    } else {
        #pragma omp for schedule(runtime)
        for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
            patch_nParticles = vecPatches( ipatch )->vecSpecies[speciesId_]->particles->size();
            vecPatches( ipatch )->vecSpecies[speciesId_]->particles->getProperty( iprop, property );
            i=0;
            j=patch_start[ipatch];
            while( i<patch_nParticles ) {
                buffer[j] = ( *property )[i];
                i++;
                j++;
            }
        }
    }
}


template<typename T>
void DiagnosticTrack::write_scalar( hid_t location, string name, T &buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist, unsigned int unit_type, unsigned int npart_global )
{
    hid_t did = H5Dcreate( location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    if( npart_global>0 ) {
        H5Dwrite( did, dtype, mem_space, file_space, transfer, &buffer );
    }
    openPMD_->writeRecordAttributes( did, unit_type );
    openPMD_->writeComponentAttributes( did, unit_type );
    H5Dclose( did );
}

template<typename T>
void DiagnosticTrack::write_component( hid_t location, string name, T &buffer, hid_t dtype, hid_t file_space, hid_t mem_space, hid_t plist, unsigned int unit_type, unsigned int npart_global )
{
    hid_t did = H5Dcreate( location, name.c_str(), dtype, file_space, H5P_DEFAULT, plist, H5P_DEFAULT );
    if( npart_global>0 ) {
        H5Dwrite( did, dtype, mem_space, file_space, transfer, &buffer );
    }
    openPMD_->writeComponentAttributes( did, unit_type );
    H5Dclose( did );
}



// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticTrack::getDiskFootPrint( int istart, int istop, Patch *patch )
{
    uint64_t footprint = 0;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    // Calculate the number of written parameters
    int nparams = 6 + nDim_particle;
    
    // Add necessary global headers approximately
    footprint += 2500;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 11250;
    
    // Add size of each parameter
    footprint += ndumps * ( uint64_t )( nparams * npart_total * 8 );
    
    return footprint;
}
