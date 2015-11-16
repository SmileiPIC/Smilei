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
#include "DiagnosticTestParticles.h"
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
    Diagnostic(Params &params, Patch* patch);
    
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
    
    //! name of the fields to dump
    std::vector<std::string> fieldsToDump;

    //! field dump output
    unsigned int fieldDump_every;
    
    //! number of time-steps for time-averaging of fields
    unsigned int ntime_step_avg;

    DiagnosticScalar scalars;
    DiagnosticProbe probes;
    DiagnosticPhaseSpace phases;

    std::vector<DiagnosticParticles*> vecDiagnosticParticles;
    std::vector<DiagnosticTestParticles*> vecDiagnosticTestParticles;
        
    
    void initScalars(Params&, Patch *patch);

    void initProbes(Params&, Patch *);
    void initPhases(Params&, Patch *);
    void initParticles(Params&);
    void initTestParticles(Params&);
    
    
    //! time-averaged field dump output
    unsigned int avgfieldDump_every;
    
    //! particle dump output
    unsigned int particleDump_every;
    
    //! scalar output every scalar_every (namelist group "DiagScalar" key "every")
    unsigned int scalar_every;
    
    double scalar_tmin;
    double scalar_tmax;
    
    
    //! list of vars for scalars to be written (empty means all)
    std::vector<std::string> scalar_vars;
    
    //! scalar output precision
    unsigned int scalar_precision;
    
    //! every for the standard pic timeloop output
    unsigned int print_every;
	
	//! vector containing phase1D structures
	std::vector<phaseStructure> vecPhase;

protected:

    std::vector<Timer> dtimer;
	
};

#endif
