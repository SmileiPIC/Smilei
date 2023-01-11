/*
 * Checkpoint.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef CHECKPOINT_H
#define CHECKPOINT_H

#include <string>
#include <vector>

#include <hdf5.h>
#include <Tools.h>

#include <H5.h>
#include "ElectroMagnBCAM_PML.h"
#include "EnvelopeBCAM_PML.h"

class Params;
class OpenPMDparams;
class SmileiMPI;
class Patch;
class SimWindow;
class ElectroMagn;
class Field;
class cField;
class Species;
class VectorPatch;
class Region;
class Collisions;

#include <csignal>

//  --------------------------------------------------------------------------------------------------------------------
//! Class Checkpoint
//  --------------------------------------------------------------------------------------------------------------------
class Checkpoint
{
public:
    Checkpoint( Params &params, SmileiMPI *smpi );
    //! Destructor for Checkpoint
    ~Checkpoint();
    
    //! Space dimension of a particle
    unsigned int nDim_particle;
    
    //! restart everything to file per processor
    void readPatchDistribution( SmileiMPI *smpi, SimWindow *simWin );
    void readRegionDistribution( Region &region );
    void restartAll( VectorPatch &vecPatches, Region &region, SmileiMPI *smpi, Params &params );
    void restartPatch( Patch *patch, Params &params, H5Read &g );
    
    //! restart field per proc
    void restartFieldsPerProc( H5Read &g, Field *field );
    void restart_cFieldsPerProc( H5Read &g, Field *field );
    
    //! load moving window parameters
    void restartMovingWindow( H5Read &f, SimWindow *simWindow );
    
    //! test before writing everything to file per processor
    //bool dump(unsigned int itime, double time, Params &params);
    void dump( VectorPatch &vecPatches, Region &region, unsigned int itime, SmileiMPI *smpi, SimWindow *simWindow, Params &params );
    // OK
    
    //! dump everything to file per processor
    void dumpAll( VectorPatch &vecPatches, Region &region, unsigned int itime,  SmileiMPI *smpi, SimWindow *simWin, Params &params );
    void dumpPatch( Patch *patch, Params &params, H5Write &g );
    
    //! incremental number of times we've done a dump
    unsigned int dump_number;
    
    //! incremental number of times we've done a dump_minutes
    unsigned int dump_minutes_times;
    
    //! this static variable is defined (in the .cpp) as false but becomes true when
    //! the signal SIGUSR1 is captured by the signal_callback_handler fnction
    static int signal_received;
    
    //! this function catches the SIGUSR1 signal and sets the signal_received to true
    static void signal_callback_handler( int signum )
    {
        MESSAGE( "----------------------------------------------" );
        MESSAGE( "Caught signal " << signum << " : dump + exit" );
        MESSAGE( "----------------------------------------------" );
        if( signum!=SIGUSR2 ) {
            signal_received = signum;
        }
    }
    
    //! start step of this run: zero if a first run, otherwise the number of the restart step
    unsigned int this_run_start_step;
    
    //! checkpoint asks to politely quit from simulation
    bool exit_asap;
    
    //! Timestep to dump everything
    unsigned int dump_step;
    
    //! Human minutes to dump everything
    double dump_minutes;
    
    //! exit once dump done
    bool exit_after_dump;
    
private:

    //! initialize the time zero of the simulation
    void initDumpCases();
    
    //! dump field per proc
    void dumpFieldsPerProc( H5Write &g, Field *field );
    void dump_cFieldsPerProc( H5Write &g, Field *field );
    
    //! dump moving window parameters
    void dumpMovingWindow( H5Write &f, SimWindow *simWindow );
    
    //! function that returns elapsed time from creator (uses private var time_reference)
    //double time_seconds();
    
    //! to dump and stop a simulation you might just check if a file named stop has been created this variable
    //! is true if since last time a file named stop appeared
    bool stop_file_seen_since_last_check;
    
    //! time of the constructor
    double time_reference;
    
    //! step at which perform a dump in case time_dump returns true
    unsigned int time_dump_step;
    
    //! keep the last keep_n_dumps dump files
    unsigned int keep_n_dumps;
    const unsigned int keep_n_dumps_max;
    
    //! write dump drectory
    std::string dump_dir;
    
    //! int deflate dump value
    int dump_deflate;
    
    std::vector<MPI_Request> dump_request;
    MPI_Status dump_status_prob;
    MPI_Status dump_status_recv;
    
    //! group checkpoint files in subdirs of file_grouping files
    unsigned int file_grouping;
    
    //! restart file
    std::string restart_file;
    
    //! dump PML in the checkpoint file 
    template <typename Tpml>
    void  dump_PML(Tpml embc, H5Write &g );
    void  dump_PML( ElectroMagnBCAM_PML *embc, H5Write &g, unsigned int imode );
    template <typename Tpml>
    void  dump_PMLenvelope(Tpml envbc, H5Write &g, unsigned int bcId );
    void  dump_PMLenvelopeAM(EnvelopeBCAM_PML *envbc, H5Write &g, unsigned int bcId );
    template <typename Tpml>
    void  restart_PML(Tpml embc, H5Read &g );
    void  restart_PML( ElectroMagnBCAM_PML *embc, H5Read &g, unsigned int imode );
    template <typename Tpml>
    void  restart_PMLenvelope(Tpml envbc, H5Read &g, unsigned int bcId );
    void  restart_PMLenvelopeAM(EnvelopeBCAM_PML *envbc, H5Read &g, unsigned int bcId );
};

#endif /* CHECKPOINT_H_ */
