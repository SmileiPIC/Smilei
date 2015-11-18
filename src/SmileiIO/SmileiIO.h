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
class Patch;
class SmileiMPI;
class Diagnostic;
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
    SmileiIO( Params& params, Diagnostic* diag, Patch* patch );
    void createFiles( Params& params, Patch* patch );
    void setFiles( hid_t masterFileId, hid_t masterFileIdAvg );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();

    
    //! Write all fields (E, B, J, rho, per species ; 10 + 4 x nspecies fields) of all time step in the same file
    void createTimeStepInSingleFileTime( int time,  Diagnostic* diag );
    void writeAllFieldsSingleFileTime( std::vector<Field*> &, int, bool );

    //! Basic Write of a field in the specified group of the global file
    virtual void writeFieldsSingleFileTime( Field* field, hid_t group_id ) = 0;

    //! Id of "Fields.h5", contains all fields per timestep
    hid_t global_file_id_;
    
    //! Id of "Fields_avg.h5", contains time-averaged fields per timestep
    hid_t global_file_id_avg;
    
    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;
    
    virtual void updatePattern( Params& params, Patch* patch ) = 0;

    //! Id of "particles-mpirank.h5", contains particles of current mpirank
    //! Disabled for now
    hid_t  partFile_id;
        
    //! Basic write field on its own file (debug)
    virtual void write( Field* field ) = 0;
    
//    template <class T> void appendTestParticles(hid_t fid, std::string name, std::vector<T> property, int nParticles, hid_t type );

//    template <class T> void appendTestParticles0( hid_t fid, std::string name, std::vector<T> property, int nParticles, hid_t type);
        
        
private:
    
    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;
    

};

#endif
