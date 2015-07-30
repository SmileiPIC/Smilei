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
    void createFiles( PicParams& params, DiagParams &diagParams, Patch* patch );
    void setFiles( hid_t masterFileId, hid_t masterFileIdAvg );
    //! Destructor for SmileiIO
    virtual ~SmileiIO();

    //! Write all fields (E, B, J, rho, per species ; 10 + 4 x nspecies fields) of all time step in the same file
    void writeAllFieldsSingleFileTime( ElectroMagn* EMfields, int itime );
    void createTimeStepInSingleFileTime( int time,  DiagParams &diagParams );

    //! Write time-averaged fields E, B) of all time step in the same file
    //! @see global_file_id_avg
    void writeAvgFieldsSingleFileTime( ElectroMagn* EMfields, int itime );
    
    //! Basic Write of a field in the specified group of the global file
    virtual void writeFieldsSingleFileTime( Field* field, hid_t group_id ) = 0;

    //! Id of "Fields.h5", contains all fields per timestep
    hid_t global_file_id_;
    
    //! Id of "Fields_avg.h5", contains time-averaged fields per timestep
    hid_t global_file_id_avg;

    //! Property list for collective dataset write, set for // IO.
    hid_t write_plist;

    //! Basic write field on its own file (debug)
    virtual void write( Field* field ) = 0;

private:
   	
	
};

#endif /* SMILEI_OUTPUT_H_ */
