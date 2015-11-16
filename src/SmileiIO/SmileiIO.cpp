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
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "ElectroMagn.h"
#include "Species.h"


using namespace std;

// static varable must be defined and initialized here
int SmileiIO::signal_received=0;

SmileiIO::SmileiIO( Params& params, Diagnostic& diag, SmileiMPI* smpi ) :
this_run_start_step(0),
dump_times(0), 
time_reference(MPI_Wtime()),
time_dump_step(0),
dump_step(0),
dump_minutes(0.0),
exit_after_dump(true),
dump_file_sequence(2),
dump_deflate(0),
restart_dir(""),
dump_request(smpi->getSize())
{
    fieldsToDump.resize(0);
    PyTools::extract("fieldsToDump", fieldsToDump);
    
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

    if (PyTools::extract("restart_dir", restart_dir) && restart_dir.back()!='/') {
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
    if (SIG_ERR == signal(SIGUSR1, SmileiIO::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR1");
    }
    if (SIG_ERR == signal(SIGUSR2, SmileiIO::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR2");
    }
    
    // ----------------------------
    // Management of global IO file
    // ----------------------------
    MPI_Info info  = MPI_INFO_NULL;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    
    write_plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(write_plist, H5FD_MPIO_INDEPENDENT);
    
    
    // Fields.h5
    // ---------
    global_file_id_  = H5Fcreate( "Fields.h5",     H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    
    // Create property list for collective dataset write: for Fields.h5
    
    H5::attr(global_file_id_, "res_time", params.res_time);
    H5::attr(global_file_id_, "every", diag.fieldDump_every);
    vector<double> my_cell_length=params.cell_length;
    my_cell_length.resize(params.nDim_field);
    H5::attr(global_file_id_, "cell_length", my_cell_length);
    H5::attr(global_file_id_, "sim_length", params.sim_length);
        
    // Fields_avg.h5
    // -------------
    global_file_id_avg = 0;
    if  (diag.ntime_step_avg!=0) {
        global_file_id_avg = H5Fcreate( "Fields_avg.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        H5::attr(global_file_id_avg, "res_time", params.res_time);
        H5::attr(global_file_id_avg, "every", diag.fieldDump_every);
        H5::attr(global_file_id_avg, "cell_length", params.cell_length);
        H5::attr(global_file_id_avg, "sim_length", params.sim_length);
    }
    
    H5Pclose(plist_id);
    
}

SmileiIO::~SmileiIO()
{
    // Management of global IO file
    H5Fclose( global_file_id_ );
    
    // Management of global IO file
    if (global_file_id_avg != 0)
        H5Fclose( global_file_id_avg );
    
    H5Pclose( write_plist );
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
//    DEBUG("[hdf] GROUP _________________________________ " << name_t.str());
    
    // Create group inside HDF5 file
    hid_t file_id;
    if( avg ) file_id = global_file_id_avg; // different file for avg fields
    else      file_id = global_file_id_;
    hid_t group_id = H5::group(file_id, name_t.str());
    
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
    
    H5Fflush( global_file_id_, H5F_SCOPE_GLOBAL );
    
}

bool SmileiIO::dump( ElectroMagn* EMfields, unsigned int itime, std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWindow, Params &params, Diagnostic &diags) {

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
        (dump_step != 0 && ((itime-this_run_start_step) % dump_step == 0)) ||
        (time_dump_step!=0 && itime==time_dump_step)) {
        dumpAll( EMfields, itime,  vecSpecies, smpi, simWindow, params, diags);
        if (exit_after_dump || ((signal_received!=0) && (signal_received != SIGUSR2))) return true;
    }
    return false;
}

void SmileiIO::dumpAll( ElectroMagn* EMfields, unsigned int itime,  std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWin, Params &params, Diagnostic &diags) {

    unsigned int num_dump=dump_times%dump_file_sequence;
    
    hid_t fid = H5Fcreate( dumpName(num_dump,smpi).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dump_times++;
    
    MESSAGEALL("Step " << itime << " : DUMP fields and particles " << dumpName(num_dump,smpi));
    
    H5::attr(fid, "Version", string(__VERSION));
    H5::attr(fid, "CommitDate", string(__COMMITDATE));
    
    H5::attr(fid, "dump_step", itime);
    
    H5::attr(fid, "Energy_time_zero", diags.scalars.Energy_time_zero);
    H5::attr(fid, "EnergyUsedForNorm", diags.scalars.EnergyUsedForNorm);
    
    dumpFieldsPerProc(fid, EMfields->Ex_);
    dumpFieldsPerProc(fid, EMfields->Ey_);
    dumpFieldsPerProc(fid, EMfields->Ez_);
    dumpFieldsPerProc(fid, EMfields->Bx_);
    dumpFieldsPerProc(fid, EMfields->By_);
    dumpFieldsPerProc(fid, EMfields->Bz_);
    if (EMfields->Ex_avg!=NULL) {
        dumpFieldsPerProc(fid, EMfields->Ex_avg);
        dumpFieldsPerProc(fid, EMfields->Ey_avg);
        dumpFieldsPerProc(fid, EMfields->Ez_avg);
        dumpFieldsPerProc(fid, EMfields->Bx_avg);
        dumpFieldsPerProc(fid, EMfields->By_avg);
        dumpFieldsPerProc(fid, EMfields->Bz_avg);
    }
    
    H5Fflush( fid, H5F_SCOPE_GLOBAL );
    
    
    unsigned int vecSpeciesSize=vecSpecies.size();
    H5::attr(fid, "species", vecSpeciesSize);
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
        ostringstream name("");
        name << setfill('0') << setw(2) << ispec;
        string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
        hid_t gid = H5::group(fid, groupName);
        
        H5::attr(gid, "partCapacity", vecSpecies[ispec]->particles.capacity());
        H5::attr(gid, "partSize", vecSpecies[ispec]->particles.size());
        
        if (vecSpecies[ispec]->particles.size()>0) {            
            for (unsigned int i=0; i<vecSpecies[ispec]->particles.Position.size(); i++) {
                ostringstream my_name("");
                my_name << "Position-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles.Position[i],dump_deflate);
            }
            
            for (unsigned int i=0; i<vecSpecies[ispec]->particles.Momentum.size(); i++) {
                ostringstream my_name("");
                my_name << "Momentum-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles.Momentum[i],dump_deflate);
            }
            
            H5::vect(gid,"Weight", vecSpecies[ispec]->particles.Weight,dump_deflate);
            H5::vect(gid,"Charge", vecSpecies[ispec]->particles.Charge,dump_deflate);
            
            if (vecSpecies[ispec]->particles.track_every) {
                H5::vect(gid,"Id", vecSpecies[ispec]->particles.Id,dump_deflate);
            }
            
            H5::vect(gid,"bmin", vecSpecies[ispec]->bmin,dump_deflate);
            H5::vect(gid,"bmax", vecSpecies[ispec]->bmax,dump_deflate);
            
        }
        H5Gclose(gid);
    }
    
    // Dump moving window status
    if (simWin!=NULL)
        H5::attr(fid, "x_moved", simWin->getXmoved());
    
    H5Fclose( fid );
    
};

void SmileiIO::dumpFieldsPerProc(hid_t fid, Field* field)
{
    H5::vect(fid, field->name, field->data_[0], field->globalDims_, H5T_NATIVE_DOUBLE, dump_deflate);
}

string SmileiIO::dumpName(unsigned int num, SmileiMPI *smpi) {
    ostringstream nameDumpTmp("");
    nameDumpTmp << "dump-" << setfill('0') << setw(1+log10(dump_file_sequence)) << num << "-" << setfill('0') << setw(1+log10(smpi->getSize())) << smpi->getRank() << ".h5" ;
    return nameDumpTmp.str();
}

void SmileiIO::restartAll( ElectroMagn* EMfields, std::vector<Species*> &vecSpecies, SmileiMPI* smpi, SimWindow* simWin, Params &params, Diagnostic &diags) {
    
    string nameDump("");

    // This will open all sequential dumps and pick the last one
    for (unsigned int num_dump=0;num_dump<dump_file_sequence; num_dump++) {
        string dump_name=restart_dir+dumpName(num_dump,smpi);
        DEBUG("Looking for restart: " << dump_name);
        ifstream f(dump_name.c_str());
        if (f.good()) {
            hid_t fid = H5Fopen( dump_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (fid >=0) {
                unsigned int itimeTmp;
                H5::getAttr(fid, "dump_step", itimeTmp);
                H5Fclose(fid);
                if (itimeTmp>this_run_start_step) {
                    this_run_start_step=itimeTmp;
                    nameDump=dump_name;
                    dump_times=num_dump;
                }
            }
        }
        f.close();
    }
    
    if (nameDump.empty()) {
        ERROR("Cannot find a valid restart file");
    }
    
    MESSAGEALL(2, " : Restarting fields and particles " << nameDump << " step=" << this_run_start_step);
    
    hid_t fid = H5Fopen( nameDump.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0) ERROR(nameDump << " is not a valid HDF5 file");
    
    string dump_version;
    H5::getAttr(fid, "Version", dump_version);
    
    string dump_date;
    H5::getAttr(fid, "CommitDate", dump_date);

    if ((dump_version != string(__VERSION)) || (dump_date != string(__COMMITDATE))) {
        WARNING ("The code version that dumped the file is " << dump_version << " of " << dump_date);
        WARNING ("                while running version is " << string(__VERSION) << " of " << string(__COMMITDATE));
    }
    
    H5::getAttr(fid, "Energy_time_zero", diags.scalars.Energy_time_zero);
    H5::getAttr(fid, "EnergyUsedForNorm", diags.scalars.EnergyUsedForNorm);
    
    
    restartFieldsPerProc(fid, EMfields->Ex_);
    restartFieldsPerProc(fid, EMfields->Ey_);
    restartFieldsPerProc(fid, EMfields->Ez_);
    restartFieldsPerProc(fid, EMfields->Bx_);
    restartFieldsPerProc(fid, EMfields->By_);
    restartFieldsPerProc(fid, EMfields->Bz_);
    
    if (EMfields->Ex_avg!=NULL) {
        restartFieldsPerProc(fid, EMfields->Ex_avg);
        restartFieldsPerProc(fid, EMfields->Ey_avg);
        restartFieldsPerProc(fid, EMfields->Ez_avg);
        restartFieldsPerProc(fid, EMfields->Bx_avg);
        restartFieldsPerProc(fid, EMfields->By_avg);
        restartFieldsPerProc(fid, EMfields->Bz_avg);
    }
    
    unsigned int vecSpeciesSize=0;
    H5::getAttr(fid, "species", vecSpeciesSize);
    if (vecSpeciesSize != vecSpecies.size()) {
        ERROR("Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<vecSpecies.size()<<")");
    }
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
    
        ostringstream name("");
        name << setfill('0') << setw(2) << ispec;
        string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
        hid_t gid = H5Gopen(fid, groupName.c_str(),H5P_DEFAULT);
        
        unsigned int partCapacity=0;
        H5::getAttr(gid, "partCapacity", partCapacity);
        vecSpecies[ispec]->particles.reserve(partCapacity,params.nDim_particle);
        
        unsigned int partSize=0;
        H5::getAttr(gid, "partSize", partSize);
        vecSpecies[ispec]->particles.initialize(partSize,params.nDim_particle);
        
        if (partSize>0) {
            for (unsigned int i=0; i<vecSpecies[ispec]->particles.Position.size(); i++) {
                ostringstream namePos("");
                namePos << "Position-" << i;
                H5::getVect(gid, namePos.str(), vecSpecies[ispec]->particles.Position[i]);
            }
            for (unsigned int i=0; i<vecSpecies[ispec]->particles.Momentum.size(); i++) {
                ostringstream namePos("");
                namePos << "Momentum-" << i;
                H5::getVect(gid, namePos.str(), vecSpecies[ispec]->particles.Momentum[i]);
            }
            H5::getVect(gid, "Weight", vecSpecies[ispec]->particles.Weight);
            H5::getVect(gid, "Charge", vecSpecies[ispec]->particles.Charge);
            if (vecSpecies[ispec]->particles.track_every) {
                H5::getVect(gid, "Id", vecSpecies[ispec]->particles.Id);
            }
            H5::getVect(gid, "bmin", vecSpecies[ispec]->bmin);
            H5::getVect(gid, "bmax", vecSpecies[ispec]->bmax);
        }
        
        H5Gclose(gid);
    }
    
    // load window status
    if (simWin!=NULL) {
        double x_moved=0.;
        H5::getAttr(fid, "x_moved", x_moved);
        simWin->setXmoved(x_moved);
    }

    H5Fclose( fid );
};

void SmileiIO::restartFieldsPerProc(hid_t fid, Field* field)
{
    hid_t did = H5Dopen (fid, field->name.c_str(),H5P_DEFAULT);
    hid_t my_sid = H5Dget_space(did);
    int datasize=H5Sget_simple_extent_ndims(my_sid);
    if (datasize != 1)
        ERROR("Data is not 1D " << datasize);

    vector<hsize_t> mydims(datasize);
    H5Sget_simple_extent_dims(my_sid,&mydims[0],NULL);
    
    if (mydims[0] != field->globalDims_)
        ERROR("Field " << field->name << " size do not match with dump " << mydims[0]  << " != " << field->globalDims_);

    H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field->data_);

    H5Sclose(my_sid);
    H5Dclose(did);
}

