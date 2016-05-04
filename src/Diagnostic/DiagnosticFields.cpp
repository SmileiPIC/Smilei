
#include <string>

#include "DiagnosticFields.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, Patch* patch, bool avg ) :
    avg(avg)
{
    string name;
    std::vector<Field*> * allFields;
    if ( ! avg ) {
        name = "DiagFields";
        filename = "new_Fields.h5";
        allFields = &(patch->EMfields->allFields);
    }
    else {
        name = "DiagFieldsAvg";
        filename = "new_Fields_avg.h5";
        allFields = &(patch->EMfields->allFields_avg);
    }
    
    if( PyTools::nComponents(name)>1 )
        ERROR("Only one "<<name<<" can be specified");
    
    std::vector<std::string> fieldsToDump(0);
    PyTools::extract("fields", fieldsToDump, name, 0);
    
    ostringstream ss("");
    fields.resize(0);
    fields_indexes.resize(0);
    bool hasfield;
    for( int i=0; i<allFields->size(); i++ ) {
        hasfield = fieldsToDump.size()==0;
        for( int j=0; j<fieldsToDump.size(); j++ ) {
            if( (*allFields)[i]->name == fieldsToDump[j] ) {
                hasfield = true;
                break;
            }
        }
        if( hasfield ) {
            ss << (*allFields)[i]->name << " ";
            fields.push_back( (*allFields)[i] );
            fields_indexes.push_back( i );
        }
    }
    MESSAGE(1,"EM fields dump "<<(avg?"(avg)":"     ")<<" :");
    MESSAGE(2, ss.str() );
    
    timeSelection = new TimeSelection( PyTools::extract_py( "every", name, 0 ), name );
    
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    if (!params.one_patch_per_MPI)
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    else
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);
    
    type_ = "Fields";
}


DiagnosticFields::DiagnosticFields( DiagnosticFields* diag, Patch* patch )
{
    avg            = diag->avg;
    filename       = diag->filename;
    fields_indexes = diag->fields_indexes;
    write_plist    = diag->write_plist;
    
    std::vector<Field*> * allFields = (avg?&(patch->EMfields->allFields_avg):&(patch->EMfields->allFields));
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
    if ( !timeSelection->theTimeIsNow(timestep) ) return false;

    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << timestep;
        
    DEBUG("[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gcreate(fileId_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);
    
    return true;
}


void DiagnosticFields::run( Patch* patch, int timestep )
{
}


void DiagnosticFields::write(int timestep)
{
    if ( !timeSelection->theTimeIsNow(timestep) ) return;
    
    // Make group name: "/0000000000", etc.
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << timestep;
    
    // Create group inside HDF5 file
    hid_t group_id = H5Gopen(fileId_, name_t.str().c_str(), H5P_DEFAULT);
    
    for (unsigned int i=0; i<fields.size(); i++) {
    
        writeFieldsSingleFileTime(fields[i], group_id );
    
        // Re-initialize average fields
        if (filename == "new_Fields_avg.h5")
            fields[i]->put_to(0.0);
    
    }
    
    H5Gclose(group_id);
    
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

}

