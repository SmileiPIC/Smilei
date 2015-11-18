/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "Params.h"
#include "Diagnostic.h"
#include "Patch.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"


using namespace std;

SmileiIO::SmileiIO( Params& params, Diagnostic *diag, Patch* patch ) : 
global_file_id_(0),
global_file_id_avg(0)
{
    fieldsToDump.resize(0);
    PyTools::extract("fieldsToDump", fieldsToDump);
    //
    // Create property list for collective dataset write.
    //
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
}


void SmileiIO::setFiles( hid_t masterFileId, hid_t masterFileIdAvg )
{
    global_file_id_ = masterFileId;
    global_file_id_avg = masterFileIdAvg;
}

void SmileiIO::createFiles( Params& params, Patch* patch)
{
    
    // ----------------------------
    // Management of global IO file
    // ----------------------------
    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    // Fields.h5
    // ---------
    global_file_id_  = H5Fcreate( "Fields.h5",     H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

    // Create property list for collective dataset write: for Fields.h5
    H5::attr(global_file_id_, "res_time", params.res_time);
    H5::attr(global_file_id_, "every", patch->Diags->fieldDump_every);

    H5::attr(global_file_id_, "res_space", params.res_space);
    vector<double> my_cell_length=params.cell_length;
    my_cell_length.resize(params.nDim_field);
    H5::attr(global_file_id_, "cell_length", my_cell_length);
    H5::attr(global_file_id_, "sim_length", params.sim_length);

    // Fields_avg.h5
    // -------------
    global_file_id_avg = 0;
    if  (patch->Diags->ntime_step_avg!=0) {
        global_file_id_avg = H5Fcreate( "Fields_avg.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        H5::attr(global_file_id_avg, "res_time", params.res_time);
        H5::attr(global_file_id_avg, "every", patch->Diags->fieldDump_every);
        H5::attr(global_file_id_avg, "res_space", params.res_space);
        H5::attr(global_file_id_avg, "cell_length", params.cell_length);
        H5::attr(global_file_id_avg, "sim_length", params.sim_length);
    }
    
    H5Pclose(plist_id);
    
}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    if (global_file_id_ != 0)
	H5Fclose( global_file_id_ );
    // Management of global IO file
    if (global_file_id_avg != 0)
        H5Fclose( global_file_id_avg );

    H5Pclose( write_plist );
}


void SmileiIO::createTimeStepInSingleFileTime( int time, Diagnostic* diag )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;
	
    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gcreate(global_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    if (diag->ntime_step_avg!=0) {
	group_id = H5Gcreate(global_file_id_avg, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group_id);
    }
 
}

// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFieldsSingleFileTime( std::vector<Field*> &fields, int time, bool avg )
{
    // Make group name: "/0000000000", etc.
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;

    // Create group inside HDF5 file
    hid_t file_id;
    if( avg ) file_id = global_file_id_avg; // different file for avg fields
    else      file_id = global_file_id_;
    //hid_t group_id = H5::group(file_id, name_t.str());
    hid_t group_id = H5Gopen(file_id, name_t.str().c_str(), H5P_DEFAULT);


    
    for (unsigned int i=0; i<fields.size(); i++) {
        if (fieldsToDump.size()==0) {
            writeFieldsSingleFileTime(fields[i], group_id );
        } else {
            for (unsigned int j=0; j<fieldsToDump.size(); j++) {
                if (fields[i]->name==fieldsToDump[j])
                    writeFieldsSingleFileTime( fields[i], group_id );
            }
        }
    }
    
    H5Gclose(group_id);
    
    H5Fflush( file_id, H5F_SCOPE_GLOBAL );
    
}

