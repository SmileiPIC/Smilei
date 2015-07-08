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

class PicParams;
class DiagParams;
class InputData;
class SmileiMPI;
class Patch;
class SimWindow;
class ElectroMagn;
class Field;
class Species;
class VectorPatch;

//  --------------------------------------------------------------------------------------------------------------------
//! Class Checkpoint
//  --------------------------------------------------------------------------------------------------------------------
class Checkpoint {
public:
    Checkpoint( PicParams& params, DiagParams &diagParams );
    //! Destructor for Checkpoint
    virtual ~Checkpoint();

    //! Space dimension of a particle
    unsigned int nDim_particle;

    //! restart everything to file per processor
    void restartAll( VectorPatch &vecPatches, unsigned int &itime,  SmileiMPI* smpi, SimWindow* simWin, PicParams &params, InputData& input_data );
    void restartPatch( ElectroMagn* EMfields,std::vector<Species*> &vecSpecies, hid_t patch_gid );

    //! restart field per proc
    void restartFieldsPerProc(hid_t fid, Field* field);

    //! load moving window parameters
    void restartMovingWindow(hid_t fid, SimWindow* simWindow);
	
    //! test before writing everything to file per processor
    bool dump(unsigned int itime, double time, PicParams &params);
    // OK
	
    //! dump everything to file per processor
    void dumpAll( VectorPatch &vecPatches, unsigned int itime,  SmileiMPI* smpi, SimWindow* simWin, PicParams &params, InputData& input_data );
    void dumpPatch( ElectroMagn* EMfields, std::vector<Species*> vecSpecies, hid_t patch_gid );

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

#endif /* CHECKPOINT_H_ */
