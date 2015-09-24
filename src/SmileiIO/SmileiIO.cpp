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
dump_times(0), 
fieldsToDump(diag.fieldsToDump),
time_reference(MPI_Wtime()),
time_dump_step(0),
dump_step(0),
dump_minutes(0.0),
exit_after_dump(true),
dump_file_sequence(2),
dump_dir(""),
dump_deflate(0),
restart_dir(""),
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

    if (PyTools::extract("dump_dir", dump_dir))
        dump_dir+="/";
    
    PyTools::extract("dump_deflate", dump_deflate);

    if (PyTools::extract("restart_dir", restart_dir))
        restart_dir+="/";
    
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

#ifdef _IO_PARTICLE
    particleSize = params.nDim_particle + 3 + 1;
    
    ostringstream name("");
    name << "particles-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
    
    // Create 1 file containing 1 dataset per Species
    partFile_id = H5Fcreate( name.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    hsize_t dims[2] = {0, particleSize};
    hsize_t max_dims[2] = {H5S_UNLIMITED, particleSize};
    hid_t file_space = H5Screate_simple(2, dims, max_dims);
    
    hid_t plist = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_layout(plist, H5D_CHUNKED);
    hsize_t chunk_dims[2] = {1, particleSize};
    H5Pset_chunk(plist, 2, chunk_dims);
    
    for (unsigned int ispec=0 ; ispec<params.params.size() ; ispec++) {
        ostringstream speciesName("");
        speciesName << params.params[ispec].type;
        
        //here we check for the presence of multiple ccurence of the same particle name... Souldn't we add a tag for each species?
        unsigned int occurrence=0;
        for (unsigned int iocc=0 ; iocc<ispec ; iocc++) {
            if (params.params[ispec].type == params.params[iocc].type)
                occurrence++;
        }
        if (occurrence>0) 
            speciesName << "_" << occurrence;
        
        hid_t did = H5Dcreate(partFile_id, speciesName.str().c_str(), H5T_NATIVE_FLOAT, file_space, H5P_DEFAULT, plist, H5P_DEFAULT);
        partDataset_id.push_back(did);
        
        hid_t tmp_space = H5Screate(H5S_SCALAR);
        
        H5::attr(partDataset_id[ispec], "Mass", params.params[ispec].mass));
        
        H5::attr(partDataset_id[ispec], "Charge",params.params[ispec].charge);
        
        H5Sclose(tmp_space);
    }
    
    H5Pclose(plist);
    H5Sclose(file_space);
    
    dims[0] = 1;
    dims[1] = particleSize;
    partMemSpace = H5Screate_simple(2, dims, NULL);
#endif
    
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
    H5::attr(global_file_id_, "res_space", params.res_space);
    H5::attr(global_file_id_, "sim_length", params.sim_length);
        
    // Fields_avg.h5
    // -------------
    global_file_id_avg = 0;
    if  (diag.ntime_step_avg!=0) {
        global_file_id_avg = H5Fcreate( "Fields_avg.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
        
        // Create property list for collective dataset write: for Fields.h5
        H5::attr(global_file_id_avg, "res_time", params.res_time);
        H5::attr(global_file_id_avg, "every", diag.fieldDump_every);
        H5::attr(global_file_id_avg, "res_space", params.res_space);
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
    
#ifdef _IO_PARTICLE
    H5Sclose(partMemSpace);
    for ( unsigned int s=0 ; s<partDataset_id.size() ; s++ )
        H5Dclose(partDataset_id[s]);
    H5Fclose(partFile_id);
#endif    
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
        if (!fieldsToDump.size()==0) {
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


// ---------------------------------------------------------------------------------------------------------------------
// Each MPI process writes is particles in its own file, data are overwritten at each call ( particles-MPI_Rank.h5 )
// In progress ...
// ---------------------------------------------------------------------------------------------------------------------
void SmileiIO::writePlasma( vector<Species*> vecSpecies, double time, SmileiMPI* smpi )
{
    
#ifdef _IO_PARTICLE
    if (smpi->isMaster()) {
        DEBUG("write species disabled");
    }
    return;
    
    for (int ispec=0 ; ispec<vecSpecies.size(); ispec++) {
        Particles* cuParticles = &(vecSpecies[ispec])->particles;
        MESSAGE(2,"write species " << ispec);
        
        for (unsigned int p=0; p<(vecSpecies[ispec])->getNbrOfParticles(); p++ ) {
            
            hid_t file_space = H5Dget_space(partDataset_id[ispec]);
            hsize_t dimsO[2];
            H5Sget_simple_extent_dims(file_space, dimsO, NULL);
            H5Sclose(file_space);
            hsize_t dims[2];
            dims[0] = dimsO[0]+1;
            dims[1] = dimsO[1];
            H5Dset_extent(partDataset_id[ispec], dims);
            
            file_space = H5Dget_space(partDataset_id[ispec]);
            hsize_t start[2];
            hsize_t count[2] = {1, particleSize};
            start[0] = dimsO[0];
            start[1] = 0;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(partDataset_id[ispec], H5T_NATIVE_DOUBLE, partMemSpace, file_space, H5P_DEFAULT, &((*cuParticles)[ p ]->position(0)));
            H5Sclose(file_space);
            
            
        } // End for p
        
    } // End for ispec
    
#endif
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
        (dump_step != 0 && (itime % dump_step == 0)) ||
        (time_dump_step!=0 && itime==time_dump_step)) {
        dumpAll( EMfields, itime,  vecSpecies, smpi, simWindow, params, diags);
        if (exit_after_dump || ((signal_received!=0) && (signal_received != SIGUSR2))) return true;
    }
    return false;
}

void SmileiIO::dumpAll( ElectroMagn* EMfields, unsigned int itime,  std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWin, Params &params, Diagnostic &diags) { 
    ostringstream nameDump("");
    nameDump << dump_dir << "dump-" << setfill('0') << setw(4) << dump_times%dump_file_sequence << "-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
    hid_t fid = H5Fcreate( nameDump.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dump_times++;
    
    MESSAGEALL("Step " << itime << " : DUMP fields and particles " << nameDump.str());    
    
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
            
            if (vecSpecies[ispec]->particles.isTestParticles) {
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
    hsize_t dims[1]={field->globalDims_};
    hid_t sid = H5Screate_simple (1, dims, NULL);    
    hid_t did = H5Dcreate (fid, field->name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0]);
    H5Dclose (did);    
    H5Sclose(sid);
}

void SmileiIO::restartAll( ElectroMagn* EMfields, unsigned int &itime,  std::vector<Species*> &vecSpecies, SmileiMPI* smpi, SimWindow* simWin, Params &params, Diagnostic &diags) {
    
    string nameDump("");

    // This will open all sequential dumps and pick the last one
    for (unsigned int i=0;i<dump_file_sequence; i++) {
        ostringstream nameDumpTmp("");
        nameDumpTmp << restart_dir << "dump-" << setfill('0') << setw(4) << i << "-" << setfill('0') << setw(4) << smpi->getRank() << ".h5" ;
        ifstream f(nameDumpTmp.str().c_str());
        if (f.good()) {
            hid_t fid = H5Fopen( nameDumpTmp.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
            if (fid >=0) {
                unsigned int itimeTmp;
                H5::getAttr(fid, "dump_step", itimeTmp);
                H5Fclose(fid);
                if (itimeTmp>itime) {
                    itime=itimeTmp;
                    nameDump=nameDumpTmp.str();
                    dump_times=i;
                }
            }
        }
        f.close();
    }
    
    if (nameDump.empty()) {
        ERROR("Cannot find a valid restart file");
    }
    
    MESSAGEALL(2, " : Restarting fields and particles " << nameDump << " step=" << itime);
    
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
        vecSpecies[ispec]->particles.initialize(partSize,params);
        
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
            if (vecSpecies[ispec]->particles.isTestParticles) {
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

