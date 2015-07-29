/*
 * Diagnostic.h
 *
 *  Created on: 3 juil. 2013
 */

#ifndef Diagnostic_H
#define Diagnostic_H

#include <vector>
#include <fstream>

#include "Interpolator.h"
#include "DiagParams.h"
#include "DiagnosticScalar.h"
#include "DiagnosticProbe.h"
#include "DiagnosticPhaseSpace.h"

class PicParams;
class SmileiMPI;
class Patch;
class VectorPatch;
class ElectroMagn;
class Species;
class SimWindow;


//! class holder for all the diagnostics: scalars, probes(0D, 1D, 2D and 3D) and phase-space
class Diagnostic {
    friend class SmileiMPI;
    friend class VectorPatch;
    friend class SimWindow;
public:
    //! creator called from main
    Diagnostic(PicParams &params,  DiagParams &dParams, Patch* patch);
    
    //! destructor
    ~Diagnostic();
    
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp);
 
    //! get a particular scalar
    double getScalar(std::string name);
    
        
protected:
    
    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace phases;
	
};

#endif
