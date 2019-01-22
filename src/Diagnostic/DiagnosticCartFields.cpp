
#include <string>

#include "DiagnosticCartFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticCartFields::DiagnosticCartFields( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &oPMD ):
    Diagnostic( oPMD )
{
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
    fn << "cartesianFields"<< ndiag <<".h5";
    filename = fn.str();
    
    // Extract the requested fields
    vector<string> fieldsToDump( 0 );
    PyTools::extract( "fields", fieldsToDump, "DiagFields", ndiag );
    
    // List all fields that are requested
    std::vector<Field *> allFields( 0 );
    ostringstream ss( "" );
    fields_indexes.resize( 0 );
    fields_names  .resize( 0 );
    hasRhoJs = false;
    // Loop fields
    for( unsigned int i=0; i<vecPatches( 0 )->EMfields->allFields.size(); i++ ) {
        string field_name = vecPatches( 0 )->EMfields->allFields[i]->name;
        bool RhoJ = field_name.at( 0 )=='J' || field_name.at( 0 )=='R';
        bool species_field = ( field_name.at( 0 )=='J' && field_name.length()>2 ) || ( field_name.at( 0 )=='R' && field_name.length()>3 );
        // If field in list of fields to dump, then add it
        if( hasField( field_name, fieldsToDump ) ) {
            ss << field_name << " ";
            fields_indexes.push_back( i );
            fields_names  .push_back( field_name );
            if( RhoJ ) {
                hasRhoJs = true;
            }
            // If field specific to a species, then allocate it
            if( species_field ) {
                Field *field = vecPatches( 0 )->EMfields->allFields[i];
                if( field->data_ != NULL ) {
                    continue;
                }
                if( field_name.substr( 0, 2 )=="Jx" ) {
                    field->allocateDims( 0, false );
                } else if( field_name.substr( 0, 2 )=="Jy" ) {
                    field->allocateDims( 1, false );
                } else if( field_name.substr( 0, 2 )=="Jz" ) {
                    field->allocateDims( 2, false );
                } else if( field_name.substr( 0, 2 )=="Rh" ) {
                    field->allocateDims();
                }
            }
        }
    }
    
    // Some output
    ostringstream p( "" );
    p << "(time average = " << time_average << ")";
    MESSAGE( 1, "Diagnostic Fields #"<<ndiag<<" "<<( time_average>1?p.str():"" )<<" :" );
    MESSAGE( 2, ss.str() );
    
    // Create new fields in each patch, for time-average storage
    vecPatches( 0 )->EMfields->allFields_avg.resize( diag_n+1 );
    if( time_average > 1 ) {
        for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ )
            vecPatches( 0 )->EMfields->allFields_avg[diag_n].push_back(
                vecPatches( 0 )->EMfields->createField( fields_names[ifield] )
            );
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
    for( unsigned int i = 0 ; i < params.nDim_field ; i++ ) {
        tot_number_of_patches /= params.global_factor[i];
    }
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( write_plist, H5FD_MPIO_COLLECTIVE );
    
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


DiagnosticCartFields::~DiagnosticCartFields()
{
    H5Pclose( write_plist );
    
    delete timeSelection;
    delete flush_timeSelection;
}


bool DiagnosticCartFields::hasField( string field_name, vector<string> fieldsToDump )
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

void DiagnosticCartFields::openFile( Params &params, SmileiMPI *smpi, bool newfile )
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

void DiagnosticCartFields::closeFile()
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



void DiagnosticCartFields::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}

bool DiagnosticCartFields::prepare( int itime )
{

    // Leave if the iteration is not the good one
    if( itime - timeSelection->previousTime( itime ) >= time_average ) {
        return false;
    }
    
    return true;
}

void DiagnosticCartFields::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers &timers )
{
    // If time-averaging, increment the average
    if( time_average>1 ) {
        #pragma omp for schedule(static)
        for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
            vecPatches( 0 )->EMfields->incrementAvgField(
                vecPatches( 0 )->EMfields->allFields[fields_indexes[ifield]], // instantaneous field
                vecPatches( 0 )->EMfields->allFields_avg[diag_n][ifield]    // averaged field
            );
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
    
    //unsigned int nPatches( vecPatches.size() );
    
    // For each field, combine all patches and write out
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
    
        // Copy the patch field to the buffer
        #pragma omp barrier
        getField( vecPatches( 0 ), ifield );
        
        #pragma omp master
        {
            // Create field dataset in HDF5
            hid_t plist_id = H5Pcreate( H5P_DATASET_CREATE );
            hid_t dset_id  = H5Dcreate( iteration_group_id, fields_names[ifield].c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT );
            H5Pclose( plist_id );
            
            // Write
            writeField( dset_id, itime );
            
            // Attributes for openPMD
            openPMD_->writeFieldAttributes( dset_id );
            openPMD_->writeRecordAttributes( dset_id, field_type[ifield] );
            openPMD_->writeFieldRecordAttributes( dset_id );
            openPMD_->writeComponentAttributes( dset_id, field_type[ifield] );
            
            // Close dataset
            H5Dclose( dset_id );
        }
    }
    
    #pragma omp master
    {
        H5Gclose( iteration_group_id );
        if( tmp_dset_id>0 ) {
            H5Dclose( tmp_dset_id );
        }
        tmp_dset_id=0;
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        }
    }
}

bool DiagnosticCartFields::needsRhoJs( int itime )
{
    return hasRhoJs && timeSelection->theTimeIsNow( itime );
}
