
#include <string>

#include "DiagnosticCartFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticCartFields::DiagnosticCartFields( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &oPMD ):
    Diagnostic( &oPMD, "DiagFields", ndiag )
{
    tmp_dset_ = NULL;
    diag_n = ndiag;
    
    filespace_firstwrite = NULL;
    memspace_firstwrite = NULL;
    filespace_reread = NULL;
    memspace_reread = NULL;
    filespace = NULL;
    memspace = NULL;
    
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
    PyTools::extractV( "fields", fieldsToDump, "DiagFields", ndiag );
    
    // List all fields that are requested
    std::vector<Field *> allFields( 0 );
    ostringstream ss( "" );
    fields_indexes.resize( 0 );
    fields_names  .resize( 0 );
    hasRhoJs = false;
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
                vecPatches( 0 )->EMfields->createField( fields_names[ifield], params )
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
    if( filespace_firstwrite ) {
        delete filespace_firstwrite;
    }
    if( memspace_firstwrite ) {
        delete memspace_firstwrite;
    }
    if( filespace_reread ) {
        delete filespace_reread;
    }
    if( memspace_reread ) {
        delete memspace_reread;
    }
    if( filespace ) {
        delete filespace;
    }
    if( memspace ) {
        delete memspace;
    }
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

void DiagnosticCartFields::openFile( Params &params, SmileiMPI *smpi )
{
    if( file_ ) {
        return;
    }
    
    // Create file
    file_ = new H5Write( filename, true );
    
    
    // Attributes for openPMD
    openPMD_->writeRootAttributes( *file_, "", "no_particles" );
    
    // Make main "data" group where everything will be stored (required by openPMD)
    data_group_ = new H5Write( file_, "data" );
    
    file_->flush();
}

void DiagnosticCartFields::closeFile()
{
    
    if( tmp_dset_ ) {
        delete tmp_dset_;
    }
    if( data_group_ ) {
        delete data_group_;
        data_group_ = NULL;
    }
    if( file_ ) {
        delete file_;
        file_ = NULL;
    }
}



void DiagnosticCartFields::init( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
{
    // create the file
    openFile( params, smpi );
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
        iteration_group_ = new H5Write( data_group_, name_t.str() );
        // Add openPMD attributes ( "basePath" )
        openPMD_->writeBasePathAttributes( *iteration_group_, itime );
        // Add openPMD attributes ( "meshesPath" )
        openPMD_->writeMeshesAttributes( *iteration_group_ );
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
            // Write
            H5Write dset = writeField( iteration_group_, fields_names[ifield], itime );
            
            // Attributes for openPMD
            openPMD_->writeFieldAttributes( dset );
            openPMD_->writeRecordAttributes( dset, field_type[ifield] );
            openPMD_->writeFieldRecordAttributes( dset );
            openPMD_->writeComponentAttributes( dset, field_type[ifield] );
            
        }
    }
    
    #pragma omp master
    {
        delete iteration_group_;
        if( tmp_dset_ ) {
            delete tmp_dset_;
        }
        tmp_dset_ = NULL;
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            file_->flush();
        }
    }
}

bool DiagnosticCartFields::needsRhoJs( int itime )
{
    return hasRhoJs && timeSelection->theTimeIsNow( itime );
}
