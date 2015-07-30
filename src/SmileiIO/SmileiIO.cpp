/*
 * SmileiIO.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "SmileiIO.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "PicParams.h"
#include "DiagParams.h"
#include "Patch.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

SmileiIO::SmileiIO( PicParams& params, DiagParams& diagParams,  Patch* patch) : 
global_file_id_(0),
global_file_id_avg(0)
{
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

void SmileiIO::createFiles( PicParams& params, DiagParams& diagParams,  Patch* patch)
{
    // Management of global IO file
    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    global_file_id_    = H5Fcreate( "Fields.h5",     H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    global_file_id_avg = 0;
    if  (diagParams.ntime_step_avg!=0)
        global_file_id_avg = H5Fcreate( "Fields_avg.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
	
    hid_t sid  = H5Screate(H5S_SCALAR);
    hid_t aid = H5Acreate (global_file_id_, "res_time", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, write_plist);
    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(params.res_time));
    H5Sclose(sid);
    H5Aclose(aid);

    sid  = H5Screate(H5S_SCALAR);
    aid = H5Acreate (global_file_id_, "every", H5T_NATIVE_UINT, sid, H5P_DEFAULT, write_plist);
    H5Awrite(aid, H5T_NATIVE_UINT, &(diagParams.fieldDump_every));
    H5Sclose(sid);
    H5Aclose(aid);
    
    hsize_t dimsPos = params.res_space.size();
    sid = H5Screate_simple(1, &dimsPos, NULL);
    aid = H5Acreate (global_file_id_, "res_space", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, write_plist);
    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(params.res_space[0]));
    H5Aclose(aid);
    H5Sclose(sid);

    dimsPos = params.sim_length.size();
    sid = H5Screate_simple(1, &dimsPos, NULL);
    vector<double> sim_length_norm=params.sim_length;
    std::transform(sim_length_norm.begin(), sim_length_norm.end(), sim_length_norm.begin(),std::bind1st(std::multiplies<double>(),1.0/params.conv_fac));

    aid = H5Acreate (global_file_id_, "sim_length", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, write_plist);
    H5Awrite(aid, H5T_NATIVE_DOUBLE, &(sim_length_norm[0]));
    H5Aclose(aid);
    H5Sclose(sid);

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


void SmileiIO::createTimeStepInSingleFileTime( int time, DiagParams &diagParams )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;
	
    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gcreate(global_file_id_, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(group_id);

    if (diagParams.ntime_step_avg!=0) {
	group_id = H5Gcreate(global_file_id_avg, name_t.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(group_id);
    }
 
}

// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int time )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;
	
    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gopen(global_file_id_, name_t.str().c_str(), H5P_DEFAULT);
	
    writeFieldsSingleFileTime( EMfields->Ex_, group_id );
    writeFieldsSingleFileTime( EMfields->Ey_, group_id );
    writeFieldsSingleFileTime( EMfields->Ez_, group_id );
    writeFieldsSingleFileTime( EMfields->Bx_m, group_id );
    writeFieldsSingleFileTime( EMfields->By_m, group_id );
    writeFieldsSingleFileTime( EMfields->Bz_m, group_id );
    writeFieldsSingleFileTime( EMfields->Jx_, group_id );
    writeFieldsSingleFileTime( EMfields->Jy_, group_id );
    writeFieldsSingleFileTime( EMfields->Jz_, group_id );
    writeFieldsSingleFileTime( EMfields->rho_, group_id );
	
    // for all species related quantities
    for (unsigned int ispec=0; ispec<EMfields->n_species; ispec++) {
        writeFieldsSingleFileTime( EMfields->rho_s[ispec], group_id );
        writeFieldsSingleFileTime( EMfields->Jx_s[ispec],  group_id );
        writeFieldsSingleFileTime( EMfields->Jy_s[ispec],  group_id );
        writeFieldsSingleFileTime( EMfields->Jz_s[ispec],  group_id );
    }
	
    H5Gclose(group_id);
	
    if (global_file_id_) H5Fflush(global_file_id_, H5F_SCOPE_GLOBAL );

}



// ---------------------------------------------------------------------------------------------------------------------
// Write all fields of all time step in the same file
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writeAvgFieldsSingleFileTime( ElectroMagn* EMfields, int time )
{
    ostringstream name_t;
    name_t.str("");
    name_t << "/" << setfill('0') << setw(10) << time;
	
    DEBUG(10,"[hdf] GROUP _________________________________ " << name_t.str());
    hid_t group_id = H5Gopen(global_file_id_avg, name_t.str().c_str(), H5P_DEFAULT);

    writeFieldsSingleFileTime( EMfields->Ex_avg, group_id );
    writeFieldsSingleFileTime( EMfields->Ey_avg, group_id );
    writeFieldsSingleFileTime( EMfields->Ez_avg, group_id );
    writeFieldsSingleFileTime( EMfields->Bx_avg, group_id );
    writeFieldsSingleFileTime( EMfields->By_avg, group_id );
    writeFieldsSingleFileTime( EMfields->Bz_avg, group_id );
	
	
    H5Gclose(group_id);

    if (global_file_id_avg) H5Fflush(global_file_id_avg, H5F_SCOPE_GLOBAL );
	
}
