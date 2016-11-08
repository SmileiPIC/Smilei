
#include <string>

#include "DiagnosticFields.h"
#include "VectorPatch.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, VectorPatch& vecPatches, int ndiag )
{
    fileId_ = 0;
    tmp_dset_id = 0;
    diag_n = ndiag;
    
    filespace_firstwrite = 0;
    memspace_firstwrite = 0;
    filespace_reread = 0;
    memspace_reread = 0;
    
    // Extract the time_average parameter
    time_average = 1;
    PyTools::extract("time_average", time_average, "DiagFields", ndiag);
    if( time_average < 1 ) time_average = 1;
    time_average_inv = 1./((double)time_average);
    
    // Define the filename
    ostringstream fn("");
    fn << "Fields"<< ndiag <<".h5";
    filename = fn.str();
    
    // Extract the requested fields
    std::vector<std::string> fieldsToDump(0);
    PyTools::extract("fields", fieldsToDump, "DiagFields", ndiag);
    
    // List all fields that are requested
    std::vector<Field*> allFields (0);
    ostringstream ss("");
    fields_indexes.resize(0);
    fields_names  .resize(0);
    bool hasfield;
    hasRhoJs = false;
    for( unsigned int i=0; i<vecPatches(0)->EMfields->allFields.size(); i++ ) {
        string field_name = vecPatches(0)->EMfields->allFields[i]->name;
        
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
        
        if( hasfield ) {
            ss << field_name << " ";
            fields_indexes.push_back( i );
            fields_names  .push_back( field_name );
            if( field_name.at(0)=='J' || field_name.at(0)=='R' )
                hasRhoJs = true;
        }
    }
    ostringstream p("");
    p << "(time average = " << time_average << ")";
    MESSAGE(1,"Diagnostic Fields #"<<ndiag<<" "<<(time_average>1?p.str():"")<<" :");
    MESSAGE(2, ss.str() );
    
    // If time average is requested, create new fields in each patch
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
        vecPatches(ipatch)->EMfields->allFields_avg.resize( diag_n+1 );
        
        if( time_average > 1 ) {
            for( unsigned int ifield=0; ifield<fields_names.size(); ifield++)
                vecPatches(ipatch)->EMfields->allFields_avg[diag_n].push_back(
                    vecPatches(ipatch)->EMfields->createField(fields_names[ifield])
                );
        }
    }
    
    // Extract the time selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagFields", ndiag ), "DiagFields" );
    
    // If the time selection contains intervals smaller than the time average, then error
    if( timeSelection->smallestInterval() < time_average )
        ERROR("Diagnostic Fields #"<<ndiag<<" has a time average too large compared to its time-selection interval ('every')");
    
    // Extract the flush time selection
    flush_timeSelection = new TimeSelection( PyTools::extract_py( "flush_every", "DiagFields", ndiag ), "DiagFields flush_every" );
    
    // Copy the total number of patches
    tot_number_of_patches = params.tot_number_of_patches;
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
}


DiagnosticFields::~DiagnosticFields()
{
    H5Pclose( write_plist );
    
    delete timeSelection;
    delete flush_timeSelection;
}


void DiagnosticFields::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if( fileId_>0 ) return;
    
    if ( newfile ) {
        // Create file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, pid);
        H5Pclose(pid);
        
        // Create property list for collective dataset write: for Fields.h5
        vector<double> my_cell_length=params.cell_length;
        my_cell_length.resize(params.nDim_field);
        H5::attr(fileId_, "res_time"    , params.res_time);
        H5::attr(fileId_, "res_space"   , params.res_space);
        H5::attr(fileId_, "cell_length" , my_cell_length);
        H5::attr(fileId_, "sim_length"  , params.sim_length);
    }
    else {
        // Open the existing file
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
    }
}

void DiagnosticFields::closeFile()
{
    if ( filespace_firstwrite>0 ) H5Sclose( filespace_firstwrite );
    if ( memspace_firstwrite >0 ) H5Sclose( memspace_firstwrite );
    if ( filespace_reread    >0 ) H5Sclose( filespace_reread );
    if ( memspace_reread     >0 ) H5Sclose( memspace_reread );
    if ( filespace           >0 ) H5Sclose( filespace );
    if ( memspace            >0 ) H5Sclose( memspace );
    if ( tmp_dset_id         >0 ) H5Dclose( tmp_dset_id );
    
    if( fileId_>0 ) H5Fclose( fileId_ );
    fileId_ = 0;
}



void DiagnosticFields::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    // create the file
    openFile( params, smpi, true );
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
}

bool DiagnosticFields::prepare( int timestep )
{
    
    // Leave if the timestep is not the good one
    if (timestep - timeSelection->previousTime(timestep) >= time_average) return false;
    
    return true;
}


void DiagnosticFields::run( SmileiMPI* smpi, VectorPatch& vecPatches, int timestep )
{
    // If time-averaging, increment the average
    if( time_average>1 ) {
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++) {
            for( unsigned int ifield=0; ifield<fields_names.size(); ifield++) {
                vecPatches(ipatch)->EMfields->incrementAvgField(
                    vecPatches(ipatch)->EMfields->allFields[fields_indexes[ifield]], // instantaneous field
                    vecPatches(ipatch)->EMfields->allFields_avg[diag_n][ifield]      // averaged field
                );
            }
        }
    }
    
    // If is not writing timestep, leave
    if (timestep - timeSelection->previousTime(timestep) != time_average-1) return;
    
    #pragma omp master
    {
        // Calculate the structure of the file depending on 1D, 2D, ...
        refHindex = (unsigned int)(vecPatches.refHindex_);
        setFileSplitting( smpi, vecPatches );
        
        // Create group for this timestep
        ostringstream name_t;
        name_t.str("");
        name_t << "/" << setfill('0') << setw(10) << timestep;
        status = H5Lexists(fileId_, name_t.str().c_str(), H5P_DEFAULT);
        if( status==0 )
           timestep_group_id = H5::group(fileId_, name_t.str().c_str());
        // Warning if file unreachable
        if( status < 0 ) WARNING("Fields diagnostics could not write");
    }
    #pragma omp barrier
    
    // Do not output diag if this timestep has already been written or if problem with file
    if( status != 0 ) return;
    
    unsigned int nPatches( vecPatches.size() );
    
    // For each field, combine all patches and write out
    for( unsigned int ifield=0; ifield < fields_indexes.size(); ifield++ ) {
        
        // Copy the patch field to the buffer
        #pragma omp barrier
        #pragma omp for schedule(static)
        for (unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++)
            getField( vecPatches(ipatch), ifield );
        
        #pragma omp master
        {
            // Create field dataset in HDF5
            hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
            hid_t dset_id  = H5Dcreate( timestep_group_id, fields_names[ifield].c_str(), H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, plist_id, H5P_DEFAULT);
            H5Pclose(plist_id);
            
            // Write
            writeField(dset_id, timestep);
            
            // Close dataset
            H5Dclose( dset_id );
        }
    }
    
    #pragma omp master
    {
        H5Gclose(timestep_group_id);
        if( tmp_dset_id>0 ) H5Dclose( tmp_dset_id );
        tmp_dset_id=0;
        if( flush_timeSelection->theTimeIsNow(timestep) ) H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
    }
}


bool DiagnosticFields::needsRhoJs(int timestep) {
    return hasRhoJs && timeSelection->theTimeIsNow(timestep);
}
