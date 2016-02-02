
#ifndef DiagnosticPhaseSpace_H
#define DiagnosticPhaseSpace_H

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

#include "DiagnosticPhase.h"

class Params;
class Patch;
class ElectroMagn;

//! mother class of all the DiagnosticPhase* it creates all the sub-diagnostics and creates and fills the hdf5 file
class DiagnosticPhaseSpace {
    friend class DiagsVectorPatch;
    friend class SimWindow;
    friend class SmileiMPI;
public:

    DiagnosticPhaseSpace(Params& params, Patch* patch);
    ~DiagnosticPhaseSpace(){};

    void run(int timestep, std::vector<Species*>& vecSpecies);
	
    void close();

    //! this vector will hold all the diagnostics created
    std::vector<DiagnosticPhase*> vecDiagPhase;

    partStruct my_part;

protected :
    //! this vector will store diagnostics to write
    std::vector<DiagnosticPhase*> vecDiagPhaseToRun;

    std::vector<DiagnosticPhase*>::iterator itDiagPhase;
	
    //! this is the hdf5 file id (we need to keep it to close at the right time)
    hid_t fileId;

};
#endif
