/*
 * SmileiIO.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef SMILEIIO_H
#define SMILEIIO_H

#include <string>
#include <vector>

#include <hdf5.h>
#include <Tools.h>

class Params;
class Diagnostic;
class SmileiMPI;
class SimWindow;
class ElectroMagn;
class Field;
class Species;

#include <csignal>

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO {
public:
    //! Create // HDF5 environment
    //! @see global_file_id_ 
    //! @see global_file_id_avg
    SmileiIO( Params& params, Diagnostic &diag, SmileiMPI* smpi );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();
    
    //! Write all fields (E, B, J, rho, per species ; 10 + 4 x nspecies fields) of all time step in the same file
    void writeAllFieldsSingleFileTime( std::vector<Field*> *, int, bool );
    
    //! Basic Write of a field in the specified group of the global file
    virtual void writeFieldsSingleFileTime( Field* field, hid_t group_id ) = 0;
    
    //! Each MPI process writes is particles in its own file
    //! Disabled for now, replaced by dump (used for restart)
    void writePlasma( std::vector<Species*> vecSpecies, double time, SmileiMPI* smpi );
    
    //! Id of "Fields.h5", contains all fields per timestep
    hid_t global_file_id_;
    
    //! Id of "Fields_avg.h5", contains time-averaged fields per timestep
    hid_t global_file_id_avg;
    
    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;
    
    //! Id of "particles-mpirank.h5", contains particles of current mpirank
    //! Disabled for now
    hid_t  partFile_id;
    
#ifdef _IO_PARTICLE
    //! Particles output in progress
    std::vector<hid_t> partDataset_id;
    //unsigned int nDatasetSpecies;
    hid_t partMemSpace;
    int particleSize;
#endif
    
    //! Basic write field on its own file (debug)
    virtual void write( Field* field ) = 0;
    
    //! restart everything to file per processor
    void restartAll( ElectroMagn* EMfields, unsigned int &itime,  std::vector<Species*> &vecSpecies, SmileiMPI* smpi, SimWindow* simWin, Params &params, Diagnostic &diags);

    //! restart field per proc
    void restartFieldsPerProc(hid_t fid, Field* field);
    
    //! load moving window parameters
    void restartMovingWindow(hid_t fid, SimWindow* simWindow);
    
    //! test before writing everything to file per processor
    bool dump(ElectroMagn* EMfields, unsigned int itime,  std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWin,  Params &params, Diagnostic &diags);

    void initWriteTestParticles(Species* species, int ispec, int itime, Params& params, SmileiMPI* smpi);
    void writeTestParticles(Species* species, int ispec, int itime, Params& params, SmileiMPI* smpi);

    template <class T> void appendTestParticles(hid_t fid, std::string name, std::vector<T> property, int nParticles, hid_t type );

    template <class T> void appendTestParticles0( hid_t fid, std::string name, std::vector<T> property, int nParticles, hid_t type);
        
    //! this static variable is defined (in the .cpp) as false but becomes true when
    //! the signal SIGUSR1 is captured by the signal_callback_handler fnction
    static int signal_received;
    
    //! this function catches the SIGUSR1 signal and sets the signal_received to true
    static void signal_callback_handler(int signum) {
        MESSAGE("----------------------------------------------");
        MESSAGE("Caught signal " << signum << " : dump + exit");
        MESSAGE("----------------------------------------------");
        if (signum!=SIGUSR2)
            signal_received = signum;
    }
    
private:
    //! incremental number of times we've done a dump
    unsigned int dump_times;
    
    //! dump everything to file per processor
    void dumpAll( ElectroMagn* EMfields, unsigned int itime,  std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWin,  Params &params, Diagnostic &diags);
    
    //! dump field per proc
    void dumpFieldsPerProc(hid_t fid, Field* field);
    
    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;
    
    //! dump moving window parameters
    void dumpMovingWindow(hid_t fid, SimWindow* simWindow);
    
    //! time of the constructor
    double time_reference;
	
    //! vector containing the stea at which perform a dump in case time_dump returns true
    unsigned int time_dump_step;
    
    //! Timestep to dump everything
    unsigned int dump_step;
    
    //! Human minutes to dump everything
    double dump_minutes;
    
    //! exit once dump done
    bool exit_after_dump;
    
    //! keep the last dump_file_sequence dump files
    unsigned int dump_file_sequence;
    
    //! write dump drectory
    std::string dump_dir;
    
    //! write dump drectory
    std::string restart_dir;
    
    std::vector<MPI_Request> dump_request;
    MPI_Status dump_status_prob;
    MPI_Status dump_status_recv;

};

#endif
