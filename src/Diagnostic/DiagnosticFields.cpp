
#include <algorithm>

#include "DiagnosticFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int ndiag, OpenPMDparams &oPMD ):
    Diagnostic( &oPMD, "DiagFields", ndiag )
{
    //MESSAGE("Starting diag field creation " );
    diag_n = ndiag;
    
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
    vector<string> allFields( vecPatches.emfields( 0 )->allFields.size() );
    for( unsigned int i=0; i<allFields.size(); i++ ) {
        allFields[i] = vecPatches.emfields( 0 )->allFields[i]->name;
    }
    // If empty list, use all fields
    if( fieldsToDump.size()==0 ) {
        fieldsToDump = allFields;
    }
    // Loop fields
    for( unsigned int j=0; j<fieldsToDump.size(); j++ ) {
        bool hasfield = false;
        for( unsigned int i=0; i<allFields.size(); i++ ) {
            // If field to dump available, then add it
            if( fieldsToDump[j] == allFields[i] ) {
                ss << allFields[i] << " ";
                fields_indexes.push_back( i );
                fields_names  .push_back( allFields[i] );
                if( allFields[i].at( 0 )=='J' || allFields[i].at( 0 )=='R' ) {
                    hasRhoJs = true;
                }
                // If field specific to a species, then allocate it
                if( params.speciesField( allFields[i] ) != "" ) {
                    vecPatches.allocateField( i, params );
                }
                hasfield = true;
                break;
            }
        }
        if( ! hasfield ) {
            ERROR_NAMELIST( 
                "Diagnostic Fields #"<<ndiag
                <<": field `"<<fieldsToDump[j]
                <<"` does not exist",
                LINK_NAMELIST + std::string("#particle-merging")
            );
        }
    }
    
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
            subgrid_stop_ .push_back( params.global_size_[isubgrid]+2 );
            subgrid_step_ .push_back( 1 );
        } else if( PyTools::py2scalar( subgrids[isubgrid], n ) ) {
            subgrid_start_.push_back( n );
            subgrid_stop_ .push_back( n + 1 );
            subgrid_step_ .push_back( 1 );
        } else if( PySlice_Check( subgrids[isubgrid] ) ) {
            Py_ssize_t start, stop, step, slicelength;
#if PY_MAJOR_VERSION == 2
            if( PySlice_GetIndicesEx( ( PySliceObject * )subgrids[isubgrid], params.global_size_[isubgrid]+1, &start, &stop, &step, &slicelength ) < 0 ) {
#else
            if( PySlice_GetIndicesEx( subgrids[isubgrid], params.global_size_[isubgrid]+1, &start, &stop, &step, &slicelength ) < 0 ) {
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
                for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
                    vecPatches( ipatch )->EMfields->allFields_avg[diag_n].push_back(
                        vecPatches( ipatch )->EMfields->createField( fields_names[ifield], params )
                    );
                }
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
    
    // Extract the datatype
    string datatype = "";
    PyTools::extract( "datatype", datatype, "DiagFields", ndiag );
    if( datatype == "double" ) {
        file_datatype_ = H5T_NATIVE_DOUBLE;
    } else if( datatype == "float" ) {
        file_datatype_ = H5T_NATIVE_FLOAT;
    } else {
        ERROR( "Diagnostic Fields #"<<ndiag<<" has an unknown datatype `"<<datatype<<"`" );
    }
    
    // Copy the total number of patches
    tot_number_of_patches = params.tot_number_of_patches;
    
    // Prepare some openPMD parameters
    field_type.resize( fields_names.size() );
    for( unsigned int ifield=0; ifield<fields_names.size(); ifield++ ) {
        string first_char = fields_names[ifield].substr( 0, 1 );
        if( first_char == "E" || first_char == "A" ) {
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
    closeFile();
    if( filespace ) {
        delete filespace;
    }
    if( memspace ) {
        delete memspace;
    }
    delete timeSelection;
    delete flush_timeSelection;
}

void DiagnosticFields::openFile( Params &, SmileiMPI *smpi )
{
    if( file_ ) {
        return;
    }
    
    // Create file
    file_ = new H5Write( filename, &smpi->world() );
    
    file_->attr( "name", diag_name_ );
    
    // Attributes for openPMD
    openPMD_->writeRootAttributes( *file_, "", "no_particles" );
    
    // Make main "data" group where everything will be stored (required by openPMD)
    data_group_ = new H5Write( file_, "data" );
    
    file_->flush();
}

void DiagnosticFields::closeFile()
{
    if( data_group_ ) {
        delete data_group_;
        data_group_ = NULL;
    }
    if( file_ ) {
        delete file_;
        file_ = NULL;
    }
}



void DiagnosticFields::init( Params &params, SmileiMPI *smpi, VectorPatch & )
{
    // create the file
    openFile( params, smpi );
}

bool DiagnosticFields::prepare( int itime )
{

    // Leave if the iteration is not the good one
    if( itime - timeSelection->previousTime( itime ) >= time_average ) {
        return false;
    }
    
    return true;
}


void DiagnosticFields::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers & )
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
        name_t << setfill( '0' ) << setw( 10 ) << itime;
        status = data_group_->has( name_t.str() );
        if( ! status ) {
            iteration_group_ = new H5Write( data_group_, name_t.str() );
            // Add openPMD attributes ( "basePath" )
            openPMD_->writeBasePathAttributes( *iteration_group_, itime );
            // Add openPMD attributes ( "meshesPath" )
            openPMD_->writeMeshesAttributes( *iteration_group_ );
        }
    }
    #pragma omp barrier
    
    // Do not output diag if this iteration has already been written or if problem with file
    if( status ) {
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
            // Write
            H5Write dset = writeField( iteration_group_, fields_names[ifield] );
            // Attributes for openPMD
            Field *f = vecPatches( 0 )->EMfields->allFields[fields_indexes[ifield]];
            vector<double> stagger( f->dims().size() );
            for( unsigned int i = 0; i < stagger.size(); i++ ) {
                stagger[i] = 0.5 * (double) f->isDual(i);
            }
            bool ends_with_m = 0 == fields_names[ifield].compare( fields_names[ifield].length()-2, 2, "_m" );
            double stagger_t = ends_with_m ? vecPatches( 0 )->EMfields->timestep*0.5 : 0.;
            openPMD_->writeFieldAttributes( dset, subgrid_start_, subgrid_step_ );
            openPMD_->writeRecordAttributes( dset, field_type[ifield], stagger_t );
            openPMD_->writeFieldRecordAttributes( dset, stagger );
            openPMD_->writeComponentAttributes( dset, field_type[ifield] );
        }
        #pragma omp barrier 
    }
    
    #pragma omp master
    {
        
        // write x_moved
        double x_moved = simWindow ? simWindow->getXmoved() : 0.;
        iteration_group_->attr( "x_moved", x_moved );
        delete iteration_group_;
        if( flush_timeSelection->theTimeIsNow( itime ) ) {
            file_->flush();
        }
    }
    #pragma omp barrier
}

bool DiagnosticFields::needsRhoJs( int itime )
{
    
    return hasRhoJs && (itime - timeSelection->previousTime( itime ) < time_average);
}

// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticFields::getDiskFootPrint( int istart, int istop, Patch * )
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
    footprint += ndumps * nfields * ( uint64_t )( total_dataset_size * (file_datatype_==H5T_NATIVE_DOUBLE?8:4) );
    
    return footprint;
}

// Calculates the intersection between a subgrid (aka slice in python) and a contiguous zone
// of the PIC grid. The zone can be a patch or a MPI patch collection.
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
            stop = std::min( zone_end, subgrid_stop );
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
            stop = std::min( zone_end, subgrid_stop );
            if( stop <= start ) {
                nsteps = 0;
            } else {
                nsteps = ( stop - start - 1 ) / subgrid_step + 1;
            }
        }
    }
}

void DiagnosticFields::findSubgridIntersection1(
    hsize_t idim,
    hsize_t &zone_offset,  // input = start of zone in full array / output = start of zone in the subgrid
    hsize_t &zone_npoints, // input = npoints of the zone in full array / output = npoints zone in the subgrid
    hsize_t &start_in_zone // output = start of subgrid in the zone
)
{
    unsigned int istart_in_zone, istart_in_file, nsteps;
    findSubgridIntersection(
        subgrid_start_[idim], subgrid_stop_[idim], subgrid_step_[idim],
        zone_offset, zone_offset + zone_npoints,
        istart_in_zone, istart_in_file, nsteps
    );
    zone_offset = istart_in_file;
    zone_npoints = nsteps;
    start_in_zone = istart_in_zone;
}
