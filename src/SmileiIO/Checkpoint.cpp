/*
 * Checkpoint.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "Checkpoint.h"

#include <sstream>
#include <iomanip>

#include <mpi.h>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "Patch.h"

using namespace std;

Checkpoint::Checkpoint( PicParams& params, DiagParams& diagParams ) : 
dump_times(0), 
dump_minutes_times(0), 
stop_file_seen_since_last_check(false)
{
    nDim_particle=params.nDim_particle;
    //particleSize = nDim_particle + 3 + 1;
    
    // 
    initDumpCases();
}


Checkpoint::~Checkpoint()
{
}


bool Checkpoint::dump( unsigned int itime, double time, PicParams &params ) { 
    if  ((params.dump_step != 0 && ((itime + 1) % params.dump_step == 0)) ||  //+1 because itime2dump receives itime + 1.
         (params.dump_minutes != 0.0 && time/60.0 > params.dump_minutes*(dump_minutes_times+1)) ||
         (params.check_stop_file && fileStopCreated())) {
        cout << "Dump global" << endl;
        if (params.dump_minutes != 0.0 && time/60.0 > params.dump_minutes*(dump_minutes_times+1) ) {
            cout << "Dump minutes" << endl;
            dump_minutes_times++;
        }	
        dump_times++;
        return true;
    }
    return false;
}

void Checkpoint::dumpAll( VectorPatch &vecPatches, unsigned int itime,  SmileiMPI* smpi, SimWindow* simWin,  PicParams &params, InputData& input_data )
{

    hid_t fid, sid, aid, tid;
		
    ostringstream nameDump("");
    nameDump << "dump-" << setfill('0') << setw(4) << (dump_times-1)%params.dump_file_sequence << "-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
    fid = H5Fcreate( nameDump.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	
    MESSAGE(2, "DUMPING fields and particles " << nameDump.str());

	
    sid  = H5Screate(H5S_SCALAR);
    tid = H5Tcopy(H5T_C_S1);
    H5Tset_size(tid, input_data.namelist.size());
    H5Tset_strpad(tid,H5T_STR_NULLTERM);
    aid = H5Acreate(fid, "Namelist", tid, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, tid, input_data.namelist.c_str());	
    H5Aclose(aid);
    H5Sclose(sid);
    H5Tclose(tid);
	
    sid = H5Screate(H5S_SCALAR);	
    aid = H5Acreate(fid, "dump_step", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_UINT, &itime);
    H5Sclose(sid);
    H5Aclose(aid);


    for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) {

	// Open a group
	ostringstream patch_name("");
	patch_name << setfill('0') << setw(6) << vecPatches(ipatch)->Hindex();
	string patchName="patch-"+patch_name.str();
	hid_t patch_gid = H5Gcreate(fid, patchName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dumpPatch( vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, patch_gid );

	// Close a group
	H5Gclose(patch_gid);

    }

    // Dump moving window status
    if (simWin!=NULL)
	dumpMovingWindow(fid, simWin);

    H5Fclose( fid );
	
}

void Checkpoint::dumpPatch( ElectroMagn* EMfields, std::vector<Species*> vecSpecies, hid_t patch_gid )
{ 
    hid_t sid, aid, gid, did;

    dumpFieldsPerProc(patch_gid, EMfields->Ex_);
    dumpFieldsPerProc(patch_gid, EMfields->Ey_);
    dumpFieldsPerProc(patch_gid, EMfields->Ez_);
    dumpFieldsPerProc(patch_gid, EMfields->Bx_);
    dumpFieldsPerProc(patch_gid, EMfields->By_);
    dumpFieldsPerProc(patch_gid, EMfields->Bz_);
    if (EMfields->Ex_avg!=NULL) {
	dumpFieldsPerProc(patch_gid, EMfields->Ex_avg);
	dumpFieldsPerProc(patch_gid, EMfields->Ey_avg);
	dumpFieldsPerProc(patch_gid, EMfields->Ez_avg);
	dumpFieldsPerProc(patch_gid, EMfields->Bx_avg);
	dumpFieldsPerProc(patch_gid, EMfields->By_avg);
	dumpFieldsPerProc(patch_gid, EMfields->Bz_avg);
    }
	
    H5Fflush( patch_gid, H5F_SCOPE_GLOBAL );
	
    sid = H5Screate(H5S_SCALAR);
    aid = H5Acreate(patch_gid, "species", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
    unsigned int vecSpeciesSize=vecSpecies.size();
    H5Awrite(aid, H5T_NATIVE_UINT, &vecSpeciesSize);
    H5Aclose(aid);
    H5Sclose(sid);

    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
	ostringstream name("");
	name << setfill('0') << setw(2) << ispec;
	string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_param.species_type;
	gid = H5Gcreate(patch_gid, groupName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				
	sid = H5Screate(H5S_SCALAR);
	aid = H5Acreate(gid, "partCapacity", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	unsigned int partCapacity=vecSpecies[ispec]->particles->capacity();
	H5Awrite(aid, H5T_NATIVE_UINT, &partCapacity);
	H5Aclose(aid);
	H5Sclose(sid);

	sid = H5Screate(H5S_SCALAR);
	aid = H5Acreate(gid, "partSize", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	unsigned int partSize=vecSpecies[ispec]->particles->size();
	H5Awrite(aid, H5T_NATIVE_UINT, &partSize);
	H5Aclose(aid);
	H5Sclose(sid);
		
	if (partSize>0) {
	    hsize_t dimsPart[1] = {vecSpecies[ispec]->getNbrOfParticles()};
			
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++) {
		ostringstream namePos("");
		namePos << "Position-" << i;
		sid = H5Screate_simple(1, dimsPart, NULL);
		did = H5Dcreate(gid, namePos.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->particles->Position[i][0]) );
		H5Dclose(did);
		H5Sclose(sid);
	    }
			
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++) {
		ostringstream namePos("");
		namePos << "Momentum-" << i;
		sid = H5Screate_simple(1, dimsPart, NULL);
		did = H5Dcreate(gid, namePos.str().c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->particles->Momentum[i][0]) );
		H5Dclose(did);
		H5Sclose(sid);
	    }
			
	    sid = H5Screate_simple(1, dimsPart, NULL);
	    did = H5Dcreate(gid, "Weight", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->particles->Weight[0]) );
	    H5Dclose(did);
	    H5Sclose(sid);
			
	    sid = H5Screate_simple(1, dimsPart, NULL);
	    did = H5Dcreate(gid, "Charge", H5T_NATIVE_SHORT, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    H5Dwrite(did, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->particles->Charge[0]) );
	    H5Dclose(did);
	    H5Sclose(sid);


	    hsize_t dimsbmin[1] = {vecSpecies[ispec]->bmin.size()};
	    sid = H5Screate_simple(1, dimsbmin, NULL);
	    did = H5Dcreate(gid, "bmin", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    H5Dwrite(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->bmin[0]) );
	    H5Dclose(did);
	    H5Sclose(sid);

	    hsize_t dimsbmax[1] = {vecSpecies[ispec]->bmax.size()};
	    sid = H5Screate_simple(1, dimsbmax, NULL);
	    did = H5Dcreate(gid, "bmax", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	    H5Dwrite(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(vecSpecies[ispec]->bmax[0]) );
	    H5Dclose(did);
	    H5Sclose(sid);
		
	    H5Gclose(gid);

	} // End if partSize

    } // End for ispec

};

void Checkpoint::dumpFieldsPerProc(hid_t fid, Field* field)
{
	hsize_t dims[1]={field->globalDims_};
	hid_t sid = H5Screate_simple (1, dims, NULL);	
	hid_t did = H5Dcreate (fid, field->name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0]);
	H5Dclose (did);	
	H5Sclose(sid);
}

void Checkpoint::dumpMovingWindow(hid_t fid, SimWindow* simWin)
{  
    double x_moved = simWin->getXmoved();

    hid_t sid = H5Screate(H5S_SCALAR);	
    hid_t aid = H5Acreate(fid, "x_moved", H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(aid, H5T_NATIVE_DOUBLE, &x_moved);
    H5Sclose(sid);
    H5Aclose(aid);

}

void Checkpoint::restartAll( VectorPatch &vecPatches, unsigned int &itime,  SmileiMPI* smpi, SimWindow* simWin, PicParams &params, InputData& input_data )
{ 
	
     string nameDump("");
	
     // This will open both dumps and pick the last one
     for (unsigned int i=0;i<params.dump_file_sequence; i++) {
	 ostringstream nameDumpTmp("");
	 nameDumpTmp << "dump-" << setfill('0') << setw(4) << i << "-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
	 ifstream f(nameDumpTmp.str().c_str());
	 if (f.good()) {
	     hid_t fid = H5Fopen( nameDumpTmp.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);			
	     hid_t aid = H5Aopen(fid, "dump_step", H5T_NATIVE_UINT);
	     unsigned int itimeTmp=0;
	     H5Aread(aid, H5T_NATIVE_UINT, &itimeTmp);	
	     H5Aclose(aid);
	     H5Fclose(fid);
	     if (itimeTmp>itime) {
		 itime=itimeTmp;
		 nameDump=nameDumpTmp.str();
		 dump_times=i;
	     }
	 }
	 f.close();
     }
	
     if (nameDump.empty()) ERROR("Cannot find a valid restart file");
	
     MESSAGE(2, "RESTARTING fields and particles " << nameDump);

     hid_t fid = H5Fopen( nameDump.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	
     hid_t aid, gid, did, sid;
	
     aid = H5Aopen(fid, "dump_step", H5T_NATIVE_UINT);
     H5Aread(aid, H5T_NATIVE_UINT, &itime);	
     H5Aclose(aid);

     for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) {

	ostringstream patch_name("");
	patch_name << setfill('0') << setw(6) << vecPatches(ipatch)->Hindex();
	string patchName="patch-"+patch_name.str();
	hid_t patch_gid = H5Gopen(fid, patchName.c_str(),H5P_DEFAULT);

	restartPatch( vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, patch_gid );

	H5Gclose(patch_gid);

     }
	
     // load window status
     if (simWin!=NULL)
	 restartMovingWindow(fid, simWin);

     H5Fclose( fid );
 };


void Checkpoint::restartPatch( ElectroMagn* EMfields,std::vector<Species*> &vecSpecies, hid_t patch_gid )
{ 
    hid_t aid, gid, did, sid;

    restartFieldsPerProc(patch_gid, EMfields->Ex_);
    restartFieldsPerProc(patch_gid, EMfields->Ey_);
    restartFieldsPerProc(patch_gid, EMfields->Ez_);
    restartFieldsPerProc(patch_gid, EMfields->Bx_);
    restartFieldsPerProc(patch_gid, EMfields->By_);
    restartFieldsPerProc(patch_gid, EMfields->Bz_);
    if (EMfields->Ex_avg!=NULL) {
	restartFieldsPerProc(patch_gid, EMfields->Ex_avg);
	restartFieldsPerProc(patch_gid, EMfields->Ey_avg);
	restartFieldsPerProc(patch_gid, EMfields->Ez_avg);
	restartFieldsPerProc(patch_gid, EMfields->Bx_avg);
	restartFieldsPerProc(patch_gid, EMfields->By_avg);
	restartFieldsPerProc(patch_gid, EMfields->Bz_avg);
    }
	
    aid = H5Aopen(patch_gid, "species", H5T_NATIVE_UINT);
    unsigned int vecSpeciesSize=0;
    H5Aread(aid, H5T_NATIVE_UINT, &vecSpeciesSize);
    H5Aclose(aid);	
    if (vecSpeciesSize != vecSpecies.size()) {
	ERROR("Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<vecSpecies.size()<<")");
    }
	
	
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
	ostringstream name("");
	name << setfill('0') << setw(2) << ispec;
	string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_param.species_type;
	gid = H5Gopen(patch_gid, groupName.c_str(),H5P_DEFAULT);
		
	aid = H5Aopen(gid, "partCapacity", H5T_NATIVE_UINT);
	unsigned int partCapacity=0;
	H5Aread(aid, H5T_NATIVE_UINT, &partCapacity);
	H5Aclose(aid);
	vecSpecies[ispec]->particles->reserve(partCapacity,nDim_particle);		

	aid = H5Aopen(gid, "partSize", H5T_NATIVE_UINT);
	unsigned int partSize=0;
	H5Aread(aid, H5T_NATIVE_UINT, &partSize);
	H5Aclose(aid);	
	vecSpecies[ispec]->particles->initialize(partSize,nDim_particle);		
		
		
	if (partSize>0) {
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++) {
		ostringstream namePos("");
		namePos << "Position-" << i;
		did = H5Dopen(gid, namePos.str().c_str(), H5P_DEFAULT);
		H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->particles->Position[i][0]);
		H5Dclose(did);
	    }
			
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++) {
		ostringstream namePos("");
		namePos << "Momentum-" << i;
		did = H5Dopen(gid, namePos.str().c_str(), H5P_DEFAULT);
		H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->particles->Momentum[i][0]);
		H5Dclose(did);
	    }
			
	    did = H5Dopen(gid, "Weight", H5P_DEFAULT);
	    H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->particles->Weight[0]);
	    H5Dclose(did);
			
	    did = H5Dopen(gid, "Charge", H5P_DEFAULT);
	    H5Dread(did, H5T_NATIVE_SHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->particles->Charge[0]);
	    H5Dclose(did);

	    did = H5Dopen(gid, "bmin", H5P_DEFAULT);
	    sid = H5Dget_space(did);
			
	    int ndims=H5Sget_simple_extent_ndims(sid);
	    vector<hsize_t> dims(ndims);
	    H5Sget_simple_extent_dims(sid,&dims[0],NULL);
			
	    vecSpecies[ispec]->bmin.resize(dims[0]);
	    H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->bmin[0]);
	    H5Dclose(did);
	    H5Sclose(sid);
			
	    did = H5Dopen(gid, "bmax", H5P_DEFAULT);
	    sid = H5Dget_space(did);
	    H5Sget_simple_extent_dims(sid,&dims[0],NULL);
			
	    vecSpecies[ispec]->bmin.resize(dims[0]);
	    H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->bmax[0]);
	    H5Dclose(did);
	    H5Sclose(sid);
	}
		
	H5Gclose(gid);
    }

}


void Checkpoint::restartFieldsPerProc(hid_t fid, Field* field)
{
	hsize_t dims[1]={field->globalDims_};
	hid_t sid = H5Screate_simple (1, dims, NULL);
	hid_t did = H5Dopen (fid, field->name.c_str(),H5P_DEFAULT);
	H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0]);
	H5Dclose (did);
	H5Sclose(sid);
}

void Checkpoint::restartMovingWindow(hid_t fid, SimWindow* simWin)
{  
    hid_t aid = H5Aopen(fid, "x_moved", H5T_NATIVE_DOUBLE);
    double x_moved=0.;
    H5Aread(aid, H5T_NATIVE_DOUBLE, &x_moved);	
    H5Aclose(aid);

    simWin->setXmoved(x_moved);

}

void Checkpoint::initDumpCases() {
	//double time_temp = MPI_Wtime();	
	//time_reference=0;
	//MPI_Allreduce(&time_temp,&time_reference,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	stop_file_seen_since_last_check=fileStopCreated();
	if (stop_file_seen_since_last_check) ERROR("File stop exists, remove it and rerun");
}

bool Checkpoint::fileStopCreated() {
	if (stop_file_seen_since_last_check) return false;
	
	int foundStopFile=0;
	ifstream f("stop");
	if (f.good()) foundStopFile=1;
	f.close();
	//int foundStopFileAll = 0;	
	//MPI_Allreduce(&foundStopFile,&foundStopFileAll,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
	
	//if (foundStopFileAll>0) {
	if (foundStopFile>0) {
		stop_file_seen_since_last_check=true;
		return true;
	} else {
		return false;
	}
}


//double Checkpoint::time_seconds() {
//	double time_temp = MPI_Wtime();	
//	double time_sec=0;
//	MPI_Allreduce(&time_temp,&time_sec,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
//	return (time_sec-time_reference);
//}





