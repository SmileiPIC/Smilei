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
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC2D_SM.h"
#include "Species.h"
#include "PatchesFactory.h"

using namespace std;

// static varable must be defined and initialized here
int Checkpoint::signal_received=0;

Checkpoint::Checkpoint( Params& params, SmileiMPI* smpi ) :
dump_number(0),
this_run_start_step(0),
exit_asap(false),
time_reference(MPI_Wtime()),
time_dump_step(0),
dump_step(0),
dump_minutes(0.0),
exit_after_dump(true),
dump_file_sequence(2),
dump_deflate(0),
restart_dir(""),
dump_request(smpi->getSize()),
file_grouping(0),
restart_number(-1)
{
    
    if( PyTools::nComponents("DumpRestart") > 0 ) {
        
        if (PyTools::extract("dump_step", dump_step, "DumpRestart")) {
            if (dump_step)
                MESSAGE(1,"Code will dump after " << dump_step << " steps");
        }
        
        if (PyTools::extract("dump_minutes", dump_minutes, "DumpRestart")) {
            if (dump_minutes>0)
                MESSAGE(1,"Code will stop after " << dump_minutes << " minutes");
        }
        
        PyTools::extract("dump_file_sequence", dump_file_sequence, "DumpRestart");
        dump_file_sequence=std::max((unsigned int)1,dump_file_sequence);
        
        PyTools::extract("exit_after_dump", exit_after_dump, "DumpRestart");
        
        PyTools::extract("dump_deflate", dump_deflate, "DumpRestart");
        
    }
    
    restart_dir = params.restart_dir;
    
    if (dump_step || dump_minutes>0) {
        if (exit_after_dump) {
            MESSAGE(1,"Code will exit after dump");
        } else {
            MESSAGE(1,"Code will continue every " << dump_step << " steps, keeping " << dump_file_sequence << " dumps");
        }
    }
    
    if (PyTools::extract("file_grouping", file_grouping, "DumpRestart") && file_grouping > 0) {
        file_grouping=std::min((unsigned int)smpi->getSize(),file_grouping);
        MESSAGE(1,"Code will group checkpoint files by "<< file_grouping << " processors");
    }

    if (PyTools::extract("restart_number", restart_number, "DumpRestart") && restart_number >= 0) {
        MESSAGE(1,"Code will restart from checkpoint number " << restart_number);
    }

    // registering signal handler
    if (SIG_ERR == signal(SIGUSR1, Checkpoint::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR1");
    }
    if (SIG_ERR == signal(SIGUSR2, Checkpoint::signal_callback_handler)) {
        WARNING("Cannot catch signal SIGUSR2");
    }
    
    nDim_particle=params.nDim_particle;
}


string Checkpoint::dumpName(unsigned int num, SmileiMPI *smpi) {
    ostringstream nameDumpTmp("");
    nameDumpTmp << "checkpoints" << PATH_SEPARATOR;
    if (file_grouping>0) {
        nameDumpTmp << setfill('0') << setw(int(1+log10(smpi->getSize()/file_grouping+1))) << smpi->getRank()/file_grouping << PATH_SEPARATOR;
    }
    
    nameDumpTmp << "dump-" << setfill('0') << setw(1+log10(dump_file_sequence)) << num << "-" << setfill('0') << setw(1+log10(smpi->getSize())) << smpi->getRank() << ".h5" ;
        return nameDumpTmp.str();
}



//bool Checkpoint::dump( unsigned int itime, double time, Params &params ) {
void Checkpoint::dump( VectorPatch &vecPatches, unsigned int itime, SmileiMPI* smpi, SimWindow* simWindow, Params &params ) {

    // check for excedeed time
    if (dump_minutes != 0.0) {
        // master checks whenever we passed the time limit
        if (smpi->isMaster() && time_dump_step==0) {
            double elapsed_time = (MPI_Wtime() - time_reference)/60.;
            if (elapsed_time > dump_minutes*(dump_number+1)) {
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
    
    if (signal_received!=0
    || (dump_step != 0 && (itime % dump_step == 0))
    || (time_dump_step!=0 && itime==time_dump_step)) {
        dumpAll( vecPatches, itime,  smpi, simWindow, params);
        if (exit_after_dump || ((signal_received!=0) && (signal_received != SIGUSR2))) {
            exit_asap=true;
        }
        signal_received=0;
        time_dump_step=0;
    }
}

void Checkpoint::dumpAll( VectorPatch &vecPatches, unsigned int itime,  SmileiMPI* smpi, SimWindow* simWin,  Params &params )
{
    
    unsigned int num_dump=dump_number%dump_file_sequence;
    
    hid_t fid = H5Fcreate( dumpName(num_dump,smpi).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    dump_number++;
    
    MESSAGEALL("Step " << itime << " : DUMP fields and particles " << dumpName(num_dump,smpi));
    
    H5::attr(fid, "Version", string(__VERSION));
    
    H5::attr(fid, "dump_step", itime);
    H5::attr(fid, "dump_number", dump_number);
    
    H5::vect( fid, "patch_count", smpi->patch_count );
    
    H5::attr(fid, "Energy_time_zero",  static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->Energy_time_zero );
    H5::attr(fid, "EnergyUsedForNorm", static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->EnergyUsedForNorm);
    H5::attr(fid, "latest_timestep",   static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->latest_timestep  );
    
    for (unsigned int ipatch=0 ; ipatch<vecPatches.size(); ipatch++) {
        
        // Open a group
        ostringstream patch_name("");
        patch_name << setfill('0') << setw(6) << vecPatches(ipatch)->Hindex();
        string patchName="patch-"+patch_name.str();
        hid_t patch_gid = H5::group(fid, patchName.c_str());
        
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
    
    dumpFieldsPerProc(patch_gid, EMfields->Ex_);
    dumpFieldsPerProc(patch_gid, EMfields->Ey_);
    dumpFieldsPerProc(patch_gid, EMfields->Ez_);
    dumpFieldsPerProc(patch_gid, EMfields->Bx_);
    dumpFieldsPerProc(patch_gid, EMfields->By_);
    dumpFieldsPerProc(patch_gid, EMfields->Bz_);
    dumpFieldsPerProc(patch_gid, EMfields->Bx_m);
    dumpFieldsPerProc(patch_gid, EMfields->By_m);
    dumpFieldsPerProc(patch_gid, EMfields->Bz_m);
    
    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name("");
        group_name << "FieldsForDiag" << idiag;
        hid_t diag_gid = H5::group(patch_gid, group_name.str());
        
        for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ )
            dumpFieldsPerProc( diag_gid, EMfields->allFields_avg[idiag][ifield] );
        
        H5Gclose(diag_gid);
    }
    
    if ( EMfields->extFields.size()>0 ) {
        for (unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if(! EMfields->emBoundCond[bcId]) continue;
            if (dynamic_cast<ElectroMagnBC1D_SM*>(EMfields->emBoundCond[bcId]) ) {
                ElectroMagnBC1D_SM* embc = static_cast<ElectroMagnBC1D_SM*>(EMfields->emBoundCond[bcId]);
                ostringstream name("");
                name << setfill('0') << setw(2) << bcId;
                string groupName="EM_boundary-species-"+name.str();
                hid_t gid = H5::group(patch_gid, groupName);
                H5::attr(gid, "Bz_xvalmin",embc->Bz_xvalmin );
                H5::attr(gid, "Bz_xvalmin",embc->Bz_xvalmax );
                H5::attr(gid, "Bz_xvalmin",embc->By_xvalmin );
                H5::attr(gid, "Bz_xvalmin",embc->By_xvalmax );
                H5Gclose(gid);
            }
            else if ( dynamic_cast<ElectroMagnBC2D_SM*>(EMfields->emBoundCond[bcId]) ) {
                ElectroMagnBC2D_SM* embc = static_cast<ElectroMagnBC2D_SM*>(EMfields->emBoundCond[bcId]);
                ostringstream name("");
                name << setfill('0') << setw(2) << bcId;
                string groupName="EM_boundary-species-"+name.str();
                hid_t gid = H5::group(patch_gid, groupName);
                H5::vect(gid, "Bx_xvalmin_Long", embc->Bx_xvalmin_Long );
                H5::vect(gid, "Bx_xvalmax_Long", embc->Bx_xvalmax_Long );
                H5::vect(gid, "By_xvalmin_Long", embc->By_xvalmin_Long );
                H5::vect(gid, "By_xvalmax_Long", embc->By_xvalmax_Long );
                H5::vect(gid, "Bz_xvalmin_Long", embc->Bz_xvalmin_Long );
                H5::vect(gid, "Bz_xvalmax_Long", embc->Bz_xvalmax_Long );
                H5::vect(gid, "Bx_yvalmin_Trans", embc->Bx_yvalmin_Trans );
                H5::vect(gid, "Bx_yvalmax_Trans", embc->Bx_yvalmax_Trans );
                H5::vect(gid, "By_yvalmin_Trans", embc->By_yvalmin_Trans );
                H5::vect(gid, "By_yvalmax_Trans", embc->By_yvalmax_Trans );
                H5::vect(gid, "Bz_yvalmin_Trans", embc->Bz_yvalmin_Trans );
                H5::vect(gid, "Bz_yvalmax_Trans", embc->Bz_yvalmax_Trans );
                H5Gclose(gid);
            }
        }
    }
    
    H5Fflush( patch_gid, H5F_SCOPE_GLOBAL );
    H5::attr(patch_gid, "species", vecSpecies.size());
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
        ostringstream name("");
        name << setfill('0') << setw(2) << ispec;
        string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
        hid_t gid = H5::group(patch_gid, groupName);
        
        H5::attr(gid, "partCapacity", vecSpecies[ispec]->particles->capacity());
        H5::attr(gid, "partSize", vecSpecies[ispec]->particles->size());
        
        if (vecSpecies[ispec]->particles->size()>0) {
            
            for (unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++) {
                ostringstream my_name("");
                my_name << "Position-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles->Position[i], dump_deflate);
            }
            
            for (unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++) {
                ostringstream my_name("");
                my_name << "Momentum-" << i;
                H5::vect(gid,my_name.str(), vecSpecies[ispec]->particles->Momentum[i], dump_deflate);
            }
            
            H5::vect(gid,"Weight", vecSpecies[ispec]->particles->Weight, dump_deflate);
            H5::vect(gid,"Charge", vecSpecies[ispec]->particles->Charge, dump_deflate);
            
            if (vecSpecies[ispec]->particles->tracked) {
                H5::vect(gid,"Id", vecSpecies[ispec]->particles->Id, dump_deflate);
            }
            
            
            H5::vect(gid,"bmin", vecSpecies[ispec]->bmin);
            H5::vect(gid,"bmax", vecSpecies[ispec]->bmax);
            
        } // End if partSize
        
        H5Gclose(gid);
        
    } // End for ispec
};


void Checkpoint::restartAll( VectorPatch &vecPatches,  SmileiMPI* smpi, SimWindow* simWin, Params &params )
{
    
    if (params.restart) {
        
        MESSAGE(1, "READING fields and particles for restart");
        
        string nameDump("");
        
        if (restart_number>=0) {
            nameDump=restart_dir+dumpName(restart_number,smpi)
        } else {
            // This will open both dumps and pick the last one
            for (unsigned int num_dump=0;num_dump<dump_file_sequence; num_dump++) {
                string dump_name=restart_dir+dumpName(num_dump,smpi);
                ifstream f(dump_name.c_str());
                if (f.good()) {
                    f.close();
                    hid_t fid = H5Fopen( dump_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
                    unsigned int stepStartTmp=0;
                    H5::getAttr(fid, "dump_step", stepStartTmp );
                    if (stepStartTmp>this_run_start_step) {
                        this_run_start_step=stepStartTmp;
                        nameDump=dump_name;
                        dump_number=num_dump;
                        H5::getAttr(fid, "dump_number", dump_number );
                    }
                    H5Fclose(fid);
                }
            }
        }
        
        if (nameDump.empty()) ERROR("Cannot find a valid restart file");
        
        MESSAGEALL(2, " : Restarting fields and particles " << nameDump << " step=" << this_run_start_step);
        
        hid_t fid = H5Fopen( nameDump.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (fid < 0) ERROR(nameDump << " is not a valid HDF5 file");
        
        string dump_version;
        H5::getAttr(fid, "Version", dump_version);
        
        string dump_date;
        H5::getAttr(fid, "CommitDate", dump_date);
        
        if (dump_version != string(__VERSION)) {
            WARNING ("The code version that dumped the file is " << dump_version);
            WARNING ("                while running version is " << string(__VERSION));
        }
        
        vector<int> patch_count(smpi->getSize());
        H5::getVect( fid, "patch_count", patch_count );
        smpi->patch_count = patch_count;
        vecPatches = PatchesFactory::createVector(params, smpi);
        
        H5::getAttr(fid, "Energy_time_zero",  static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->Energy_time_zero );
        H5::getAttr(fid, "EnergyUsedForNorm", static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->EnergyUsedForNorm);
        H5::getAttr(fid, "latest_timestep",   static_cast<DiagnosticScalar*>(vecPatches.globalDiags[0])->latest_timestep  );
        
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
    }
    
}


void Checkpoint::restartPatch( ElectroMagn* EMfields,std::vector<Species*> &vecSpecies, Params& params, hid_t patch_gid )
{
    restartFieldsPerProc(patch_gid, EMfields->Ex_);
    restartFieldsPerProc(patch_gid, EMfields->Ey_);
    restartFieldsPerProc(patch_gid, EMfields->Ez_);
    restartFieldsPerProc(patch_gid, EMfields->Bx_);
    restartFieldsPerProc(patch_gid, EMfields->By_);
    restartFieldsPerProc(patch_gid, EMfields->Bz_);
    restartFieldsPerProc(patch_gid, EMfields->Bx_m);
    restartFieldsPerProc(patch_gid, EMfields->By_m);
    restartFieldsPerProc(patch_gid, EMfields->Bz_m);
    
    // Fields required for DiagFields
    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++ ) {
        ostringstream group_name("");
        group_name << "FieldsForDiag" << idiag;
        hid_t diag_gid = H5Gopen(patch_gid, group_name.str().c_str(),H5P_DEFAULT);
        
        for( unsigned int ifield=0; ifield<EMfields->allFields_avg[idiag].size(); ifield++ )
            restartFieldsPerProc( diag_gid, EMfields->allFields_avg[idiag][ifield] );
        
        H5Gclose(diag_gid);
    }
    
    if ( EMfields->extFields.size()>0 ) {
        for (unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
            if(! EMfields->emBoundCond[bcId]) continue;
            if (dynamic_cast<ElectroMagnBC1D_SM*>(EMfields->emBoundCond[bcId]) ) {
                ElectroMagnBC1D_SM* embc = static_cast<ElectroMagnBC1D_SM*>(EMfields->emBoundCond[bcId]);
                ostringstream name("");
                name << setfill('0') << setw(2) << bcId;
                string groupName="EM_boundary-species-"+name.str();
                hid_t gid = H5Gopen(patch_gid, groupName.c_str(),H5P_DEFAULT);
                H5::getAttr(gid, "Bz_xvalmin", embc->Bz_xvalmin );
                H5::getAttr(gid, "Bz_xvalmin", embc->Bz_xvalmax );
                H5::getAttr(gid, "Bz_xvalmin", embc->By_xvalmin );
                H5::getAttr(gid, "Bz_xvalmin", embc->By_xvalmax );
                H5Gclose(gid);
                
            }
            else if ( dynamic_cast<ElectroMagnBC2D_SM*>(EMfields->emBoundCond[bcId]) ) {
                ElectroMagnBC2D_SM* embc = static_cast<ElectroMagnBC2D_SM*>(EMfields->emBoundCond[bcId]);
                ostringstream name("");
                name << setfill('0') << setw(2) << bcId;
                string groupName="EM_boundary-species-"+name.str();
                hid_t gid = H5Gopen(patch_gid, groupName.c_str(),H5P_DEFAULT);
                H5::getVect(gid, "Bx_xvalmin_Long", embc->Bx_xvalmin_Long );
                H5::getVect(gid, "Bx_xvalmax_Long", embc->Bx_xvalmax_Long );
                H5::getVect(gid, "By_xvalmin_Long", embc->By_xvalmin_Long );
                H5::getVect(gid, "By_xvalmax_Long", embc->By_xvalmax_Long );
                H5::getVect(gid, "Bz_xvalmin_Long", embc->Bz_xvalmin_Long );
                H5::getVect(gid, "Bz_xvalmax_Long", embc->Bz_xvalmax_Long );
                H5::getVect(gid, "Bx_yvalmin_Trans", embc->Bx_yvalmin_Trans );
                H5::getVect(gid, "Bx_yvalmax_Trans", embc->Bx_yvalmax_Trans );
                H5::getVect(gid, "By_yvalmin_Trans", embc->By_yvalmin_Trans );
                H5::getVect(gid, "By_yvalmax_Trans", embc->By_yvalmax_Trans );
                H5::getVect(gid, "Bz_yvalmin_Trans", embc->Bz_yvalmin_Trans );
                H5::getVect(gid, "Bz_yvalmax_Trans", embc->Bz_yvalmax_Trans );
                H5Gclose(gid);
            }
        }
    }
    
    unsigned int vecSpeciesSize=0;
    H5::getAttr(patch_gid, "species", vecSpeciesSize );
    
    if (vecSpeciesSize != vecSpecies.size()) {
        ERROR("Number of species differs between dump (" << vecSpeciesSize << ") and namelist ("<<vecSpecies.size()<<")");
    }
    
    
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
        ostringstream name("");
        name << setfill('0') << setw(2) << ispec;
        string groupName="species-"+name.str()+"-"+vecSpecies[ispec]->species_type;
        hid_t gid = H5Gopen(patch_gid, groupName.c_str(),H5P_DEFAULT);
        
        unsigned int partCapacity=0;
        H5::getAttr(gid, "partCapacity", partCapacity );
        vecSpecies[ispec]->particles->reserve(partCapacity,nDim_particle);
        
        unsigned int partSize=0;
        H5::getAttr(gid, "partSize", partSize );
        vecSpecies[ispec]->particles->initialize(partSize,nDim_particle);
        
        
        if (partSize>0) {
            for (unsigned int i=0; i<vecSpecies[ispec]->particles->Position.size(); i++) {
                ostringstream namePos("");
                namePos << "Position-" << i;
                H5::getVect(gid,namePos.str(),vecSpecies[ispec]->particles->Position[i]);
            }
            
            for (unsigned int i=0; i<vecSpecies[ispec]->particles->Momentum.size(); i++) {
                ostringstream namePos("");
                namePos << "Momentum-" << i;
                H5::getVect(gid,namePos.str(),vecSpecies[ispec]->particles->Momentum[i]);
            }
            
            H5::getVect(gid,"Weight",vecSpecies[ispec]->particles->Weight);
            
            H5::getVect(gid,"Charge",vecSpecies[ispec]->particles->Charge);
            
            H5::getVect(gid,"Id",vecSpecies[ispec]->particles->Id);

            H5::getVect(gid,"bmin",vecSpecies[ispec]->bmin,true);
            H5::getVect(gid,"bmax",vecSpecies[ispec]->bmax,true);

        }
        
        H5Gclose(gid);
    }
    
}

void Checkpoint::dumpFieldsPerProc(hid_t fid, Field* field)
{
    hsize_t dims[1]={field->globalDims_};
    hid_t sid = H5Screate_simple (1, dims, NULL);
    hid_t did = H5Dcreate (fid, field->name.c_str(), H5T_NATIVE_DOUBLE, sid, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &field->data_[0]);
    H5Dclose (did);
    H5Sclose(sid);
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

void Checkpoint::dumpMovingWindow(hid_t fid, SimWindow* simWin)
{
    H5::attr(fid, "x_moved", simWin->getXmoved());
    H5::attr(fid, "n_moved", simWin->getNmoved());
    
}
void Checkpoint::restartMovingWindow(hid_t fid, SimWindow* simWin)
{

    double x_moved=0.;
    H5::getAttr(fid, "x_moved", x_moved );
    simWin->setXmoved(x_moved);
    
    unsigned int n_moved=0;
    H5::getAttr(fid, "n_moved", n_moved );
    simWin->setNmoved(n_moved);
    
}






