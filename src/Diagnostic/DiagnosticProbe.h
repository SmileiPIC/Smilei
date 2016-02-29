/*
-----------------------------------------------------------------------
PROBE DIAGNOSTICS                   - Mickael & Julien ?? - 2014
                                    - Modified by F Perez - 04/2015
-----------------------------------------------------------------------
See doc for help

*/


#ifndef DiagnosticProbe_H
#define DiagnosticProbe_H

#include <cmath>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <hdf5.h>

#include "Tools.h"

#include "Species.h"
#include "Interpolator.h"
#include "Particles.h"
#include "TimeSelection.h"

class Params;
class Patch;
class ElectroMagn;
class Field2D;

//! this class holds the point probe
class DiagnosticProbe {
    friend class SmileiMPI;
 
public:
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe(Params &params, Patch* patch);
    //DiagnosticProbe();

    void createFile();
    void setFile(hid_t masterFileId, Patch* patch, Params& params);
    void waitSetFile(Params& params);

    void setFile(hid_t masterFileId);
    //! function to close the file
    void close();

    void writePositionIn( Params &params );
    void writePositions(int probe_id, int ndim_Particles, int probeDim, hid_t group_id );
    
    ~DiagnosticProbe();
    
    //! run all probes
    void run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    
    void compute(int probe_id, unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    void write(int probe_id, unsigned int timestep, hid_t group_id);

    //! return name of the probe based on its number
    std::string probeName(int p);

    //! vector containing the time selections at which to calculate each probe
    std::vector<TimeSelection*> timeSelection;

    //! hdf5 file ID
    hid_t fileId;

    // rank of the cpu (from smpi) -> patch->hindex
    const unsigned int cpuRank;
    
    double dt;
    
    //! fake particles acting as probes
    std::vector<Particles> probeParticles;
    
    //! each probe will write in a buffer
    std::vector< Field2D* > probesArray;
    
    std::vector<int> probesStart;

    //! number of fake particles for each probe diagnostic
    std::vector<unsigned int> nPart_total;
    int nProbeTot;
    
    //! Number of fields to save
    std::vector<int> nFields;
    
    //! List of fields to save
    std::vector<std::vector<std::string>> fieldname;
    
    //! Indices in the output array where each field goes
    std::vector<std::vector<unsigned int>> fieldlocation;
    
protected:
    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    //! J local fields for the projector
    LocalFields Jloc_fields;
    //! Rho local field for the projector
    double Rloc_fields;
    
    //! memory size of a probe should be 6 = Exyz + Bxyz
    const int probeSize;   

private:
    std::vector<MPI_Request> rsend;
    std::vector<MPI_Request> rrecv;

};
#endif
