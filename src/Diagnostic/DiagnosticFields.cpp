
#include <string>

#include "DiagnosticFields.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, Patch* patch, int ndiag )
{
    // Extract the time_average parameter
    time_average = 1;
    PyTools::extract("time_average", time_average, "DiagFields", ndiag);
    if( time_average < 1 )
        time_average = 1;
    
    // Verify that only one diag of this type exists
    int tavg;
    for( int idiag=0; idiag<ndiag; idiag++ ) {
        tavg = 1;
        PyTools::extract("time_average", tavg, "DiagFields", idiag);
        if( tavg < 1 )
            tavg = 1;
        if( tavg*time_average == 1 || (tavg>1 && time_average>1) )
            ERROR("Cannot have two DiagFields with time_average "<<(tavg==1?"=":">")<<" 1");
    }
    
    // Define the filename and get the list of fields
    std::vector<Field*> * allFields;
    if ( time_average==1 ) {
        filename = "Fields.h5";
        allFields = &(patch->EMfields->allFields);
    }
    else {
        filename = "Fields_avg.h5";
        allFields = &(patch->EMfields->allFields_avg);
    }
    
    // Extract the requested fields
    std::vector<std::string> fieldsToDump(0);
    PyTools::extract("fields", fieldsToDump, "DiagFields", ndiag);
    
    // List all fields that are requested
    ostringstream ss("");
    fields.resize(0);
    fields_indexes.resize(0);
    bool hasfield;
    for( int i=0; i<allFields->size(); i++ ) {
        string field_name = (*allFields)[i]->name;
        if( field_name.find("_avg") < string::npos ) field_name.erase(field_name.find("_avg"));
        
        if( fieldsToDump.size()==0 ) {
            hasfield = true;
        } else {
            hasfield = false;
            for( int j=0; j<fieldsToDump.size(); j++ ) {
                if( field_name == fieldsToDump[j] ) {
                    hasfield = true;
                    break;
                }
            }
        }
        
        if( hasfield ) {
            ss << field_name << " ";
            fields.push_back( (*allFields)[i] );
            fields_indexes.push_back( i );
        }
    }
    MESSAGE(1,"EM fields dump "<<(time_average>1?"(avg)":"     ")<<" :");
    MESSAGE(2, ss.str() );
    
    // Extract the time selection
    timeSelection = new TimeSelection( PyTools::extract_py( "every", "DiagFields", ndiag ), "DiagFields" );
    
    // Prepare the property list for HDF5 output
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    if (!params.one_patch_per_MPI)
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    else
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
    type_ = "Fields";
}


DiagnosticFields::DiagnosticFields( DiagnosticFields* diag, Patch* patch )
{
    time_average   = diag->time_average;
    filename       = diag->filename;
    fields_indexes = diag->fields_indexes;
    write_plist    = diag->write_plist;
    
    std::vector<Field*> * allFields = (time_average>1?&(patch->EMfields->allFields_avg):&(patch->EMfields->allFields));
    fields.resize(0);
    for( int i=0; i<fields_indexes.size(); i++ )
        fields.push_back( (*allFields)[fields_indexes[i]] );
    
    timeSelection = new TimeSelection( diag->timeSelection );
    
    type_ = "Fields";
}


DiagnosticFields::~DiagnosticFields()
{
    // Management of global IO file
    if (fileId_ != 0)
        H5Fclose(fileId_ );

    H5Pclose( write_plist );
    
    delete timeSelection;
}


void DiagnosticFields::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if ( newfile ) {
        // ----------------------------
        // Management of global IO file
        // ----------------------------
        MPI_Info info  = MPI_INFO_NULL;
        hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
        
        // Fields.h5
        // ---------
        fileId_  = H5Fcreate( filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        vector<double> my_cell_length=params.cell_length;
        my_cell_length.resize(params.nDim_field);
        H5::attr(fileId_, "res_time", params.res_time);
        H5::attr(fileId_, "res_space", params.res_space);
        H5::attr(fileId_, "cell_length", my_cell_length);
        H5::attr(fileId_, "sim_length", params.sim_length);
        
        H5Pclose(plist_id);
    }
    else {
        hid_t pid = H5Pcreate(H5P_FILE_ACCESS);
        H5Pset_fapl_mpio(pid, MPI_COMM_WORLD, MPI_INFO_NULL);
        fileId_ = H5Fopen( filename.c_str(), H5F_ACC_RDWR, pid );
        H5Pclose(pid);
    }
    
}


void DiagnosticFields::closeFile()
{
    H5Fclose(fileId_);
}


bool DiagnosticFields::prepare( int timestep )
{
    // Get the previous selected time
    int previousTime = timeSelection->previousTime(timestep);
    
    // Leave if the timestep is not the good one
    if (timestep - previousTime >= time_average) return false;
    
    // Prepare HDF5 group if at right timestep
    if (timestep - previousTime == time_average-1) {
        ostringstream name_t;
        name_t.str("");
        name_t << "/" << setfill('0') << setw(10) << timestep;
            
        DEBUG("[hdf] GROUP _________________________________ " << name_t.str());
        hid_t group_id = H5Gcreate(fileId_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);
    }
    
    return true;
}


void DiagnosticFields::run( Patch* patch, int timestep )
{
    // If time-averaging, increment the average
    if( time_average>1 )
        patch->EMfields->incrementAvgFields(timestep);
    
    // If writing timestep, then write
    if (timestep - timeSelection->previousTime(timestep) == time_average-1) {
        // Make group name: "/0000000000", etc.
        ostringstream name_t;
        name_t.str("");
        name_t << "/" << setfill('0') << setw(10) << timestep;
        
        // Create group inside HDF5 file
        hid_t group_id = H5Gopen(fileId_, name_t.str().c_str(), H5P_DEFAULT);
        
        for (unsigned int i=0; i<fields.size(); i++) {
            
            writeField(fields[i], group_id );
            
            // Re-initialize average fields
            if (time_average>1) fields[i]->put_to(0.0);
        }
        
        H5Gclose(group_id);
        
        H5Fflush( fileId_, H5F_SCOPE_GLOBAL );
        
        patch->EMfields->restartRhoJs();
    }
}


void DiagnosticFields::write(int timestep)
{
}

