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

class PicParams;
class SmileiMPI;
class Patch;
class DiagParams;
class ElectroMagn;
class Field2D;

//! this class holds the point probe
class DiagnosticProbe {
    
public:
    friend class SmileiMPI;
    
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe(PicParams &params, DiagParams &diagParams, Patch* patch);

    void createFile(DiagParams &diagParams);
    void setFile(hid_t masterFileId, Patch* patch, PicParams& params, DiagParams &diagParams);
    void setFile(hid_t masterFileId);

    void writePositionIn( PicParams &params, DiagParams &diagParams );
    void writePositions(int probe_id, int ndim_Particles, int probeDim, hid_t group_id );
    
    ~DiagnosticProbe();
    
    //! run all probes
    void run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    
    void compute(int probe_id, unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    void write(int probe_id, unsigned int timestep, hid_t group_id);

    //! return name of the probe based on its number
    std::string probeName(int p);

    //! function to close the file
    void close();

    //! vector containing the timesteps at which calculate each probe
    std::vector<unsigned int> every;
    //! hdf5 file ID
    hid_t fileId;

protected:
    
    // rank of the cpu (from smpi) -> patch->hindex
    const unsigned int cpuRank;
    
    //! fake particles acting as probes
    std::vector<Particles> probeParticles;
    
    //! each probe will write in a buffer
    std::vector< Field2D* > probesArray;
    std::vector< int > probesStart;
    int nProbeTot;
    //int nDim;
    
    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    
    //! J local fields for the projector
    LocalFields Jloc_fields;
    
    //! memory size of a probe should be 6 = Exyz + Bxyz
    const int probeSize;
    

 private:
};
#endif
