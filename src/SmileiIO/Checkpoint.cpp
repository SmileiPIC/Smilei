/*
 * Checkpoint.cpp
 *
 *  Created on: 3 juil. 2013
 */

#include "Checkpoint.h"

#include <sstream>
#include <iomanip>
#include <string>

#include <mpi.h>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "VectorPatch.h"
#include "Diagnostic.h"

using namespace std;

// static varable must be defined and initialized here
int Checkpoint::signal_received=0;

Checkpoint::Checkpoint( Params& params, SmileiMPI* smpi ) : 
this_run_start_step(0),
dump_times(0), 
time_reference(MPI_Wtime()),
time_dump_step(0),
dump_step(0),
dump_minutes(0.0),
exit_after_dump(true),
dump_file_sequence(2),
dump_request(smpi->getSize())
{
    if (PyTools::extract("dump_step", dump_step)) {
        if (dump_step)
            MESSAGE(1,"Code will dump after " << dump_step << " steps");
    }
    
    if (PyTools::extract("dump_minutes", dump_minutes)) {
        if (dump_minutes>0)
            MESSAGE(1,"Code will stop after " << dump_minutes << " minutes");
    }
    
    PyTools::extract("dump_file_sequence", dump_file_sequence);
    dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);
    
    PyTools::extract("exit_after_dump", exit_after_dump);

    PyTools::extract("dump_deflate", dump_deflate);

    if (PyTools::extract("restart_dir", restart_dir) && restart_dir.at(restart_dir.length()-1)!='/') {
        restart_dir+="/";
    }

    if (dump_step || dump_minutes>0) {
        if (exit_after_dump) {
            MESSAGE(1,"Code will exit after dump");
        } else {
            MESSAGE(1,"Code will continue every " << dump_step << " steps, keeping " << dump_file_sequence << " dumps");
        }
    }
    
    // registering signal handler
    if (SIG_ERR == signal(SIGUSR1, Checkpoint::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR1");
    }
    if (SIG_ERR == signal(SIGUSR2, Checkpoint::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR2");
    }
/*    // one of these below should be the soft linit signal for loadlever
    if (SIG_ERR == signal(SIGXCPU, SmileiIO::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGXCPU");
    }
    if (SIG_ERR == signal(SIGTERM, SmileiIO::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGTERM");
    }
    if (SIG_ERR == signal(SIGINT, SmileiIO::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGINT");
    }
*/

    nDim_particle=params.nDim_particle;
    //particleSize = nDim_particle + 3 + 1;
    
    // 
}


Checkpoint::~Checkpoint()
{
}


//bool Checkpoint::dump( unsigned int itime, double time, Params &params ) { 
bool Checkpoint::dump( VectorPatch &vecPatches, unsigned int itime, SmileiMPI* smpi, SimWindow* simWindow, Params &params, Diagnostic* diags ) { 

    // check for excedeed time 
    if (dump_minutes != 0.0) {
        // master checks whenever we passed the time limit
        if (smpi->isMaster() && time_dump_step==0) {
            double elapsed_time = (MPI_Wtime() - time_reference)/60.;
            if (elapsed_time > dump_minutes*(dump_times+1)) {
                time_dump_step = itime+1; // we will dump at next timestep (in case non-master already passed)
                MESSAGE("Reached time limit : " << elapsed_time << " minutes. Dump timestep : " << time_dump_step );
                // master does a non-blocking send
                for (unsigned int dest=0; dest < (unsigned int) smpi->getSize(); dest++) {
                    MPI_Isend(&time_dump_step,1,MPI_UNSIGNED,dest,SMILEI_COMM_DUMP_TIME,smpi->SMILEI_COMM_WORLD,&dump_request[dest]);
                }
            }
        } else { // non master nodes receive the time_dump_step (non-blocking)
            int todump=0;
            MPI_Iprobe(0,SMILEI_COMM_DUMP_TIME,MPI_COMM_WORLD,&todump,&dump_status_prob);
            if (todump) {
                MPI_Recv(&time_dump_step,1,MPI_UNSIGNED,0,SMILEI_COMM_DUMP_TIME,smpi->SMILEI_COMM_WORLD,&dump_status_recv);
            }
        }
    }
    
    if (signal_received!=0 ||
        (dump_step != 0 && (itime % dump_step == 0)) ||
        (time_dump_step!=0 && itime==time_dump_step)) {
        dumpAll( vecPatches, itime,  smpi, simWindow, params, diags);
        if (exit_after_dump || ((signal_received!=0) && (signal_received != SIGUSR2))) return true;
    }
    return false;
    

    /*if  ((params.dump_step != 0 && ((itime + 1) % params.dump_step == 0)) ||  //+1 because itime2dump receives itime + 1.
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
    return false;*/
}

void Checkpoint::dumpAll( VectorPatch &vecPatches, unsigned int itime,  SmileiMPI* smpi, SimWindow* simWin,  Params &params, Diagnostic* diags )
{

    hid_t fid, sid, aid, tid;
		
    ostringstream nameDump("");
    nameDump << "dump-" << setfill('0') << setw(4) << dump_times%dump_file_sequence << "-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
    fid = H5Fcreate( nameDump.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dump_times++;

    MESSAGEALL("Step " << itime << " : DUMP fields and particles " << nameDump.str());    

    H5::attr(fid, "Version", string(__VERSION));
    H5::attr(fid, "CommitDate", string(__COMMITDATE));
    
    H5::attr(fid, "dump_step", itime);
    
    H5::attr(fid, "Energy_time_zero", diags->scalars.Energy_time_zero);
    H5::attr(fid, "EnergyUsedForNorm", diags->scalars.EnergyUsedForNorm);

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
	
    H5::attr(patch_gid, "species", vecSpecies.size());    

    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
	ostringstream name("");
	name << setfill('0') << setw(2) << ispec;
	string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
	hid_t gid = H5::group(patch_gid, groupName);
		
	sid = H5Screate(H5S_SCALAR);
	aid = H5Acreate(gid, "partCapacity", H5T_NATIVE_UINT, sid, H5P_DEFAULT, H5P_DEFAULT);
	unsigned int partCapacity=vecSpecies[ispec]->particles->capacity();
	H5Awrite(aid, H5T_NATIVE_UINT, &partCapacity);
	H5Aclose(aid);
	H5Sclose(sid);

        H5::attr(gid, "partCapacity", vecSpecies[ispec]->particles->capacity());
        H5::attr(gid, "partSize", vecSpecies[ispec]->particles->size());
		
	if (vecSpecies[ispec]->particles->size()>0) {
	    hsize_t dimsPart[1] = {vecSpecies[ispec]->getNbrOfParticles()};
			
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++) {
		ostringstream my_name("");
                my_name << "Position-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles->Position[i]);
	    }
			
	    for (unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++) {
		ostringstream my_name("");
                my_name << "Momentum-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles->Momentum[i]);
	    }
			
            H5::vect(gid,"Weight", vecSpecies[ispec]->particles->Weight);
            H5::vect(gid,"Charge", vecSpecies[ispec]->particles->Charge);

            if (vecSpecies[ispec]->particles->track_every) {
                H5::vect(gid,"Id", vecSpecies[ispec]->particles->Id);
            }


            H5::vect(gid,"bmin", vecSpecies[ispec]->bmin);
            H5::vect(gid,"bmax", vecSpecies[ispec]->bmax);
		
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
    H5::attr(fid, "x_moved", simWin->getXmoved());
    H5::attr(fid, "n_moved", simWin->getNmoved());

}

void Checkpoint::restartAll( VectorPatch &vecPatches, unsigned int &itime,  SmileiMPI* smpi, SimWindow* simWin, Params &params )
{ 
	
     string nameDump("");
	
     // This will open both dumps and pick the last one
     for (unsigned int i=0;i<dump_file_sequence; i++) {
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
	
     MESSAGEALL(2, " : Restarting fields and particles " << nameDump << " step=" << this_run_start_step);

     hid_t fid = H5Fopen( nameDump.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
     if (fid < 0) ERROR(nameDump << " is not a valid HDF5 file");
    
     string dump_version;
     H5::getAttr(fid, "Version", dump_version);
    
     string dump_date;
     H5::getAttr(fid, "CommitDate", dump_date);

     if ((dump_version != string(__VERSION)) || (dump_date != string(__COMMITDATE))) {
         WARNING ("The code version that dumped the file is " << dump_version << " of " << dump_date);
         WARNING ("                while running version is " << string(__VERSION) << " of " << string(__COMMITDATE));
     }

     H5::getAttr(fid, "Energy_time_zero", vecPatches.Diags->scalars.Energy_time_zero);
     H5::getAttr(fid, "EnergyUsedForNorm", vecPatches.Diags->scalars.EnergyUsedForNorm);
	
     hid_t aid, gid, did, sid;
	
     aid = H5Aopen(fid, "dump_step", H5T_NATIVE_UINT);
     H5Aread(aid, H5T_NATIVE_UINT, &itime);	
     H5Aclose(aid);



     for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) {

	ostringstream patch_name("");
	patch_name << setfill('0') << setw(6) << vecPatches(ipatch)->Hindex();
	string patchName="patch-"+patch_name.str();
	hid_t patch_gid = H5Gopen(fid, patchName.c_str(),H5P_DEFAULT);

	restartPatch( vecPatches(ipatch)->EMfields, vecPatches(ipatch)->vecSpecies, params, patch_gid );

	H5Gclose(patch_gid);

     }
	
     // load window status
     if (simWin!=NULL)
	 restartMovingWindow(fid, simWin);

     H5Fclose( fid );
 };


void Checkpoint::restartPatch( ElectroMagn* EMfields,std::vector<Species*> &vecSpecies, Params& params, hid_t patch_gid )
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
	string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
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
	    
            if (vecSpecies[ispec]->particles->track_every) {
                did = H5Dopen(gid, "Id", H5P_DEFAULT);
                H5Dread(did, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vecSpecies[ispec]->particles->Id[0]);
                H5Dclose(did);
            }

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

    aid = H5Aopen(fid, "n_moved", H5T_NATIVE_UINT);
    unsigned int n_moved=0;
    H5Aread(aid, H5T_NATIVE_UINT, &n_moved);	
    H5Aclose(aid);
    
    simWin->setXmoved(x_moved);
    simWin->setNmoved(n_moved);

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





