
#include <string>

#include "DiagnosticFields.h"

using namespace std;

DiagnosticFields::DiagnosticFields( Params &params, SmileiMPI* smpi, Patch* patch, int diagId )
{
    string name, name_every, time_name;
    if (diagId == 0) { // diagId = 0, std
        name = "fieldsToDump";
        name_every = "fieldDump_every";
        time_name = "Fields";
        filename = "new_Fields.h5";
        fields_ = patch->EMfields->allFields;
    }
    else if (diagId == 1) { // diagId = 1, avg
        name = "avgfieldsToDump";
        name_every = "avgfieldDump_every";
        time_name = "Average fields";
        filename = "Fields_avg.h5";
        fields_ = patch->EMfields->allFields_avg;
    }
    else 
        ERROR( "This should not exist !!!" );

    fieldsToDump.resize(0);
    PyTools::extract(name, fieldsToDump);

    //
    // Create property list for collective dataset write.
    //
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    if (!params.simu_is_cartesian)
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    else
        H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_COLLECTIVE);

    timeSelection    = new TimeSelection( PyTools::extract_py( name_every ), time_name );

    type_ = "Fields";
}


DiagnosticFields::DiagnosticFields()
{
}


DiagnosticFields::~DiagnosticFields()
{
    // Management of global IO file
    if (fileId_ != 0)
        H5Fclose(fileId_ );

    H5Pclose( write_plist );
}


void DiagnosticFields::openFile( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, bool newfile )
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


bool DiagnosticFields::prepare( Patch* patch, int timestep )
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

    for (unsigned int i=0; i<fields_.size(); i++) {

        if (fieldsToDump.size()==0) { // Write all fields
            writeFieldsSingleFileTime(fields_[i], group_id );
        }
        else { // Write selected fields only
            for (unsigned int j=0; j<fieldsToDump.size(); j++) {
                if (fields_[i]->name==fieldsToDump[j]) {
                    writeFieldsSingleFileTime( fields_[i], group_id );
                }
            }
        }

        // Re-initialize average fields
        if (filename == "Fields_avg.h5")
            fields_[i]->put_to(0.0);

    } // END for i
    
    H5Gclose(group_id);
    
    H5Fflush( fileId_, H5F_SCOPE_GLOBAL );

}

