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

class PicParams;
class DiagParams;
class InputData;
class SmileiMPI;
class Patch;
class SimWindow;
class ElectroMagn;
class Field;
class Species;

//  --------------------------------------------------------------------------------------------------------------------
//! Class SmileiIO
//  --------------------------------------------------------------------------------------------------------------------
class SmileiIO {
public:
    //! Create // HDF5 environment
    //! @see global_file_id_ 
    //! @see global_file_id_avg
    SmileiIO( PicParams& params, DiagParams &diagParams, Patch* patch );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();

    //! Write all fields (E, B, J, rho, per species ; 10 + 4 x nspecies fields) of all time step in the same file
    void writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int itime );
    
    //! Write time-averaged fields E, B) of all time step in the same file
    //! @see global_file_id_avg
    void writeAvgFieldsSingleFileTime( ElectroMagn* EMfields, int itime );
    
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

    //! Space dimension of a particle
    unsigned int nDim_particle;

    //! Basic write field on its own file (debug)
    virtual void write( Field* field ) = 0;

    //! restart everything to file per processor
    void restartAll( ElectroMagn* EMfields, unsigned int &itime,  std::vector<Species*> &vecSpecies, SmileiMPI* smpi, SimWindow* simWin, PicParams &params, InputData& input_data);

    //! restart field per proc
    void restartFieldsPerProc(hid_t fid, Field* field);

    //! load moving window parameters
    void restartMovingWindow(hid_t fid, SimWindow* simWindow);
	
    //! test before writing everything to file per processor
    bool dump(ElectroMagn* EMfields, unsigned int itime, double time,  std::vector<Species*> vecSpecies, SimWindow* simWin,  PicParams &params, InputData& input_data);
	
    //! dump everything to file per processor
    void dumpAll( ElectroMagn* EMfields, unsigned int itime,  std::vector<Species*> vecSpecies, SmileiMPI* smpi, SimWindow* simWin,  PicParams &params, InputData& input_data);

    //! incremental number of times we've done a dump
    unsigned int dump_times;

    //! incremental number of times we've done a dump_minutes
    unsigned int dump_minutes_times;

private:
    
    //! initialize the time zero of the simulation 
    void initDumpCases();
	
    //! dump field per proc
    void dumpFieldsPerProc(hid_t fid, Field* field);

    //! dump moving window parameters
    void dumpMovingWindow(hid_t fid, SimWindow* simWindow);

    //! time of the constructor
    //double time_reference;
	
    //! function that returns elapsed time from creator (uses private var time_reference)
    //double time_seconds();
	
    //! to dump and stop a simulation you might just check if a file named stop has been created this variable
    //! is true if since last time a file named stop appeared
    bool stop_file_seen_since_last_check;

    //! function that checks if file named "stop" exists;
    bool fileStopCreated();
	
	
};

#endif /* SMILEI_OUTPUT_H_ */
