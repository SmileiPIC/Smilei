#ifndef DiagnosticProbe1D_H
#define DiagnosticProbe1D_H

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
class DiagParams;
class ElectroMagn;

//! this class holds the point probe
class DiagnosticProbe1D {

public:

    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe1D(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticProbe1D();

    void run(int timestep, unsigned int numDiag, ElectroMagn* EMfields, Interpolator* interp);

    //! return name of the probe based on its number
    std::string probeName(int p);
    
    //! function to close the file
    void close();

    //! vector containing the timesteps at which calculate each probe1D
    std::vector<unsigned int> every;

private:
    //! copy of the smpi pointer
    SmileiMPI* smpi_;

    //! fake particles acting as probes
    std::vector<Particles> probeParticles;
    
    //! each probe will be associated with a proc
    std::vector<int> probeId;
    
    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    
    //! memory size of a probe should be 7 = timestep + Exyz + Bxyz
    const int probeSize;

    //! hdf5 file ID
    hid_t fileId;
    
    //! hdf5 group ID of each probe
    std::vector<hid_t> groupId;

};
#endif
