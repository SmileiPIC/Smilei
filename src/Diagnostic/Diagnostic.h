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
#include "DiagnosticScalar.h"
#include "DiagnosticProbe.h"
#include "DiagnosticPhaseSpace.h"
#include "DiagnosticParticles.h"
#include "DiagnosticTrackParticles.h"
#include "Timer.h"

class Params;
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
    Diagnostic(Params &params, Patch* patch, SmileiMPI* smpi);
    
    //! destructor
    ~Diagnostic(){};
    
    //! print timers
    void printTimers(Patch *patch, double tottime);
    
    //! close all diags
    void closeAll(Patch* patch);

    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp);
 
    //! get a particular scalar
    double getScalar(std::string name);

    std::vector<Timer> dtimer;

    DiagnosticScalar scalars;
    DiagnosticProbe probes;
    DiagnosticPhaseSpace phases;

    std::vector<DiagnosticParticles*> vecDiagnosticParticles;
    std::vector<DiagnosticTrackParticles*> vecDiagnosticTrackParticles;
    
    //! Time selection for field dumps
    TimeSelection * field_timeSelection;
    
    //! Time selection for average field dumps
    TimeSelection * avgfield_timeSelection;
    
    //! number of time-steps for time-averaging of fields
    unsigned int ntime_step_avg;
    
    //! every for the standard pic timeloop output
    unsigned int print_every;
	
};

#endif
