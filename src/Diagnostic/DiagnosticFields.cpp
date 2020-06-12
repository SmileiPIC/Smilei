
#include <string>

#include "DiagnosticFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &oPMD ):
    Diagnostic( &oPMD, "DiagFields", ndiag )
{
    //MESSAGE("Starting diag field creation " );
    fileId_ = 0;
    data_group_id = 0;
    tmp_dset_id = 0;
    diag_n = ndiag;
    
    filespace_firstwrite = 0;
    memspace_firstwrite = 0;
    filespace_reread = 0;
    memspace_reread = 0;
    
    // Extract the time_average parameter
    time_average = 1;
    PyTools::extract( "time_average", time_average, "DiagFields", ndiag );
    if( time_average < 1 ) {
        time_average = 1;
    }
    time_average_inv = 1./( ( double )time_average );
    
    // Define the filename
    ostringstream fn( "" );
    fn << "Fields"<< ndiag <<".h5";
    filename = fn.str();
    
    // Extract the requested fields
    vector<string> fieldsToDump( 0 );
    PyTools::extractV( "fields", fieldsToDump, "DiagFields", ndiag );
   
    //Avoid modes repetition in the namelist by interpreting field quantity as all modes of this quantity
    if (params.geometry == "AMcylindrical") {
        vector<string> fieldsToAdd( 0 );
        for( unsigned int ifield = 0; ifield < fieldsToDump.size(); ifield++ ){
            if( fieldsToDump[ifield].find("_mode_") ==  std::string::npos ) {
                for( unsigned int imode = 1; imode < params.nmodes; imode ++ ){
                    fieldsToAdd.push_back(fieldsToDump[ifield]+"_mode_"+to_string(imode));
                }
                fieldsToDump[ifield] = fieldsToDump[ifield] + "_mode_0" ;
            }
        }
        for( unsigned int ifield = 0; ifield < fieldsToAdd.size(); ifield++ ){
            fieldsToDump.push_back(fieldsToAdd[ifield]);
        }
    }
 
    // List all fields that are requested
    ostringstream ss( "" );
    fields_indexes.resize( 0 );
    fields_names  .resize( 0 );
    hasRhoJs = false;
    //MESSAGE("HNA");
    // Loop fields
    for( unsigned int i=0; i<vecPatches.emfields( 0 )->allFields.size(); i++ ) {
        string field_name = vecPatches.emfields( 0 )->allFields[i]->name;
        
        // If field in list of fields to dump, then add it
        if( hasField( field_name, fieldsToDump ) ) {
            ss << field_name << " ";
            fields_indexes.push_back( i );
            fields_names  .push_back( field_name );
            if( field_name.at( 0 )=='J' || field_name.at( 0 )=='R' ) {
                hasRhoJs = true;
            }
            // If field specific to a species, then allocate it
            if( params.speciesField( field_name ) != "" ) {
                vecPatches.allocateField( i, params );
            }
        }
    }
    //MESSAGE("possible bug");
    
    // Extract subgrid info
    PyObject *subgrid = PyTools::extract_py( "subgrid", "DiagFields", ndiag );
    // Make a vector of all subgrids
    vector<PyObject *> subgrids;
    if( subgrid == Py_None ) {
        subgrids.resize( params.nDim_field, Py_None );
    } else if( ! PySequence_Check( subgrid ) ) {
        subgrids.push_back( subgrid );
    } else {
        Py_ssize_t ns = PySequence_Length( subgrid );
        for( Py_ssize_t is=0; is<ns; is++ ) {
            subgrids.push_back( PySequence_Fast_GET_ITEM( subgrid, is ) );
        }
    }
    Py_DECREF( subgrid );
    // Verify the number of subgrids
    unsigned int nsubgrid = subgrids.size();
    if( nsubgrid != params.nDim_field ) {
        ERROR( "Diagnostic Fields #"<<ndiag<<" `subgrid` containing "<<nsubgrid<<" axes whereas simulation dimension is "<<params.nDim_field );
    }
    // Check each subgrid is a slice, and save the slice boundaries
    for( unsigned int isubgrid=0; isubgrid<nsubgrid; isubgrid++ ) {
        unsigned int n;
        if( subgrids[isubgrid] == Py_None ) {
            subgrid_start_.push_back( 0 );
            subgrid_stop_ .push_back( params.n_space_global[isubgrid]+2 );
            subgrid_step_ .push_back( 1 );
        } else if( PyTools::py2scalar( subgrids[isubgrid], n ) ) {
            subgrid_start_.push_back( n );
            subgrid_stop_ .push_back( n + 1 );
            subgrid_step_ .push_back( 1 );
        } else if( PySlice_Check( subgrids[isubgrid] ) ) {
            Py_ssize_t start, stop, step, slicelength;
#if PY_MAJOR_VERSION == 2
            if( PySlice_GetIndicesEx( ( PySliceObject * )subgrids[isubgrid], params.n_space_global[isubgrid]+1, &start, &stop, &step, &slicelength ) < 0 ) {
#else
            if( PySlice_GetIndicesEx( subgrids[isubgrid], params.n_space_global[isubgrid]+1, &start, &stop, &step, &slicelength ) < 0 ) {
#endif
                PyTools::checkPyError();
                ERROR( "Diagnostic Fields #"<<ndiag<<" `subgrid` axis #"<<isubgrid<<" not understood" );
            }
            subgrid_start_.push_back( start );
            subgrid_stop_ .push_back( stop );
            subgrid_step_ .push_back( step );
            if( slicelength < 1 ) {
                ERROR( "Diagnostic Fields #"<<ndiag<<" `subgrid` axis #"<<isubgrid<<" is an empty selection" );
            }
        } else {
            ERROR( "Diagnostic Fields #"<<ndiag<<" `subgrid` axis #"<<isubgrid<<" must be an integer or a slice" );
        }
    }
    
    // Some output
    ostringstream p( "" );
    p << "(time average = " << time_average << ")";
    MESSAGE( 1, "Diagnostic Fields #"<<ndiag<<" "<<( time_average>1?p.str():"" )<<" :" );
    MESSAGE( 2, ss.str() );
    
    // Create new fields in each patch, for time-average storage
    if( ! smpi->test_mode ) {
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
            vecPatches( ipatch )->EMfields->allFields_avg.resize( diag_n+1 );
            if( time_average > 1 ) {
                for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ )
                    vecPatches( ipatch )->EMfields->allFields_avg[diag_n].push_back(
                        vecPatches( ipatch )->EMfields->createField( fields_names[ifield],params )
                    );
            }
        }
    }
    
    // Extract the time selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagFields", ndiag ), "DiagFields" );
    
    // If the time selection contains intervals smaller than the time average, then error
    if( timeSelection->smallestInterval() < time_average ) {
        ERROR( "Diagnostic Fields #"<<ndiag<<" has a time average too large compared to its time-selection interval ('every')" );
    }
    
    // Extract the flush time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py( "flush_every", "DiagFields", ndiag ), "DiagFields flush_every" );
    
    // Copy the total number of patches
    tot_number_of_patches = params.tot_number_of_patches;
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( write_plist, H5FD_MPIO_COLLECTIVE );
    dcreate = H5Pcreate( H5P_DATASET_CREATE );
    
    // Prepare some openPMD parameters
    field_type.resize( fields_names.size() );
    for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
        string first_char = fields_names[ifield].substr( 0, 1 );
        if( first_char == "E" ) {
            field_type[ifield] = SMILEI_UNIT_EFIELD;
        } else if( first_char == "B" ) {
            field_type[ifield] = SMILEI_UNIT_BFIELD;
        } else if( first_char == "J" ) {
            field_type[ifield] = SMILEI_UNIT_CURRENT;
        } else if( first_char == "R" ) {
            field_type[ifield] = SMILEI_UNIT_DENSITY;
        } else {
            ERROR( " impossible field name " );
        }
    }
}


DiagnosticFields::~DiagnosticFields()
{
    H5Pclose( write_plist );
    H5Pclose( dcreate );
    
    delete timeSelection;
    delete flush_timeSelection;
}


bool DiagnosticFields::hasField( string field_name, vector<string> fieldsToDump )
{
    bool hasfield;
    if( fieldsToDump.size()==0 ) {
        hasfield = true;
    } else {
        hasfield = false;
        for( unsigned int j=0; j<fieldsToDump.size(); j++ ) {
            if( field_name == fieldsToDump[j] ) {
                hasfield = true;
                break;
            }
        }
    }
    return hasfield;
}

void DiagnosticFields::openFile( Params &params, SmileiMPI *smpi, bool newfile )
{
    if( fileId_>0 ) {
        return;
    }
    
    if( newfile ) {
        // Create file
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid );
        H5Pclose( pid );
        
        H5::attr( fileId_, "name", diag_name_ );
        
        // Attributes for openPMD
        openPMD_->writeRootAttributes( fileId_, "", "no_particles" );
        
        // Make main "data" group where everything will be stored (required by openPMD)
        data_group_id = H5::group( fileId_, "data" );
    } else {
        // Open the existing file
        hid_t pid = H5Pcreate( H5P_FILE_ACCESS );
        H5Pset_fapl_mpio( pid, MPI_COMM_WORLD, MPI_INFO_NULL );
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose( pid );
        data_group_id = H5Gopen( fileId_, "data", H5P_DEFAULT );
    }
}

void DiagnosticFields::closeFile()
{
    if( filespace_firstwrite>0 ) {
        H5Sclose( filespace_firstwrite );
    }
    if( memspace_firstwrite >0 ) {
        H5Sclose( memspace_firstwrite );
    }
    if( filespace_reread    >0 ) {
        H5Sclose( filespace_reread );
    }
    if( memspace_reread     >0 ) {
        H5Sclose( memspace_reread );
    }
    if( filespace           >0 ) {
        H5Sclose( filespace );
    }
    if( memspace            >0 ) {
        H5Sclose( memspace );
    }
    if( tmp_dset_id         >0 ) {
        H5Dclose( tmp_dset_id );
    }
    
    if( data_group_id>0 ) {
        H5Gclose( data_group_id );
    }
    data_group_id = 0;
    if( fileId_>0 ) {
        H5Fclose( fileId_ );
    }
    fileId_ = 0;
}



void DiagnosticFields::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}

bool DiagnosticFields::prepare( int itime )
{

    // Leave if the iteration is not the good one
    if( itime - timeSelection->previousTime( itime ) >= time_average ) {
        return false;
    }
    
    return true;
}


void DiagnosticFields::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers )
{
    // If time-averaging, increment the average
    if( time_average>1 ) {
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
            for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
                vecPatches( ipatch )->EMfields->incrementAvgField(
                    vecPatches( ipatch )->EMfields->allFields[fields_indexes[ifield]], // instantaneous field
                    vecPatches( ipatch )->EMfields->allFields_avg[diag_n][ifield]    // averaged field
                );
            }
        }
    }
    
    // If is not writing iteration, leave
    if( itime - timeSelection->previousTime( itime ) != time_average-1 ) {
        return;
    }
    
    #pragma omp master
    {
        // Calculate the structure of the file depending on 1D, 2D, ...
        refHindex = ( unsigned int )( vecPatches.refHindex_ );
        setFileSplitting( smpi, vecPatches );
        
        // Create group for this iteration
        ostringstream name_t;
        name_t.str( "" );
        name_t << setfill( '0' ) << setw( 10 ) << itime;
        status = H5Lexists( data_group_id, name_t.str().c_str(), H5P_DEFAULT );
        if( status==0 ) {
            iteration_group_id = H5::group( data_group_id, name_t.str().c_str() );
        }
        // Warning if file unreachable
        if( status < 0 ) {
            WARNING( "Fields diagnostics could not write" );
        }
        // Add openPMD attributes ( "basePath" )
        openPMD_->writeBasePathAttributes( iteration_group_id, itime );
        // Add openPMD attributes ( "meshesPath" )
        openPMD_->writeMeshesAttributes( iteration_group_id );
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( status != 0 ) {
        return;
    }
    
    unsigned int nPatches( vecPatches.size() );
    
    // For each field, combine all patches and write out
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
    
        // Copy the patch field to the buffer
        #pragma omp barrier
        #pragma omp for schedule(static)
        for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
            getField( vecPatches( ipatch ), ifield );
        }
        
        #pragma omp master
        {
            // Create field dataset in HDF5
            hid_t dset_id  = H5Dcreate( iteration_group_id, fields_names[ifield].c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, dcreate, H5P_DEFAULT );
            
            // Write
            writeField( dset_id, itime );
            
            // Attributes for openPMD
            openPMD_->writeFieldAttributes( dset_id, subgrid_start_, subgrid_step_ );
            openPMD_->writeRecordAttributes( dset_id, field_type[ifield] );
            openPMD_->writeFieldRecordAttributes( dset_id );
            openPMD_->writeComponentAttributes( dset_id, field_type[ifield] );
            
            // Close dataset
            H5Dclose( dset_id );
        }
        #pragma omp barrier

    }
    
    #pragma omp master
    {
        // write x_moved
        double x_moved = simWindow ? simWindow->getXmoved() : 0.;
        H5::attr( iteration_group_id, "x_moved", x_moved );
        
        H5Gclose( iteration_group_id );
        if( tmp_dset_id>0 ) {
            H5Dclose( tmp_dset_id );
        }
        tmp_dset_id=0;
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        }
    }
    #pragma omp barrier
}

bool DiagnosticFields::needsRhoJs( int itime )
{
    return hasRhoJs && timeSelection->theTimeIsNow( itime );
}

// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticFields::getDiskFootPrint( int istart, int istop, Patch *patch )
{
    uint64_t footprint = 0;
    uint64_t nfields = fields_indexes.size();
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    // Add necessary global headers approximately
    footprint += 2500;
    
    // Add necessary timestep headers approximately
    footprint += ndumps * 2200;
    
    // Add necessary field headers approximately
    footprint += ndumps * nfields * 1200;
    
    // Add size of each field
    footprint += ndumps * nfields * ( uint64_t )( total_dataset_size * 8 );
    
    return footprint;
}

// Calculates the intersection between a subgrid (aka slice in python) and a contiguous zone
// of the PIC grid. The zone can be a patch or a MPI region.
void DiagnosticFields::findSubgridIntersection(
    unsigned int subgrid_start,
    unsigned int subgrid_stop,
    unsigned int subgrid_step,
    unsigned int zone_begin,
    unsigned int zone_end,
    unsigned int &istart_in_zone,  // index since the zone start that is its first intersection with subgrid
    unsigned int &istart_in_file,  // index of this first intersection in the file (or equivalently in subgrid)
    unsigned int &nsteps  // Number of intersecting elements
)
{
    unsigned int start, stop;
    // If the zone begins before the subgrid
    if( zone_begin <= subgrid_start ) {
        istart_in_zone = subgrid_start - zone_begin;
        istart_in_file = 0;
        if( zone_end <= subgrid_start ) {
            nsteps = 0;
        } else {
            stop = min( zone_end, subgrid_stop );
            if( stop <= subgrid_start ) {
                stop = subgrid_start + 1;
            }
            nsteps = ( stop - subgrid_start - 1 ) / subgrid_step + 1;
        }
    } else {
        if( zone_begin >= subgrid_stop ) {
            istart_in_zone = 0;
            istart_in_file = 0;
            nsteps = 0;
        } else {
            istart_in_file = ( zone_begin - subgrid_start - 1 ) / subgrid_step + 1;
            start = subgrid_start + istart_in_file * subgrid_step;
            istart_in_zone = start - zone_begin;
            stop = min( zone_end, subgrid_stop );
            if( stop <= start ) {
                nsteps = 0;
            } else {
                nsteps = ( stop - start - 1 ) / subgrid_step + 1;
            }
        }
    }
}
