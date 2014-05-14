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
class DiagParams;
class ElectroMagn;

//! this class holds the point probe
class DiagnosticProbe {
    
public:
    
    //! the creator need both sim parameters params and the diagnostic parameter diagParams
    DiagnosticProbe(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi, unsigned int dim_probe):
    dimProbe(dim_probe+2), smpi_(smpi), probeSize(6), fileId(0) {open();};
    
    ~DiagnosticProbe();
    
    //! run all probes
    void runAll(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    
    //! run a single probe
    void run(unsigned int numDiag, ElectroMagn* EMfields, Interpolator* interp);
    
    //! return name of the probe based on its number
    std::string probeName(int p);

    //! function to open the file
    void open();

    //! prpare the hdf5 dataset and attributes
    void addProbe(unsigned int np, std::vector<double> partPos, std::vector<unsigned int> vecNumber);

    //! function to close the file
    void close();

    //! vector containing the timesteps at which calculate each probe
    std::vector<unsigned int> every;

protected:
    //!dimension of the probe in the hdf5file (probe0D -> 2(time,fields), probe1D -> 3(time,fields,number_probes_1st_direction) etc...)
    unsigned int dimProbe;
    
    //! copy of the smpi pointer
    SmileiMPI* smpi_;
    
    //! fake particles acting as probes
    std::vector<Particles> probeParticles;
    
    //! each probe will be associated with a proc
    std::vector<std::vector<int> >probeId;
    
    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    
    //! memory size of a probe should be 6 = Exyz + Bxyz
    int probeSize;
    
    //! hdf5 file ID
    hid_t fileId;
    
};
#endif
