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
class SmileiMPI;
class ElectroMagn;
class Species;


//! class holder for all the diagnostics: scalars, probes(0D, 1D, 2D and 3D) and phase-space
class Diagnostic {

public:
    //! creator called from main
    Diagnostic(Params&, std::vector<Species*>& vecSpecies, SmileiMPI *smpi);
    
    //! destructor
    ~Diagnostic(){};
    
    //! print timers
    void printTimers(SmileiMPI *smpi, double tottime);
    
    //! close all diags
    void closeAll(SmileiMPI* smpi);

    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp, SmileiMPI *smpi);
 
    //! get a particular scalar
    double getScalar(std::string name);

    std::vector<Timer> dtimer;
    
    DiagnosticScalar scalars;
    DiagnosticProbe probes;
	DiagnosticPhaseSpace phases;
    std::vector<DiagnosticParticles*> vecDiagnosticParticles;
    std::vector<DiagnosticTestParticles*> vecDiagnosticTestParticles;
        
    void initScalars(Params&, SmileiMPI *smpi);
    void initProbes(Params&, SmileiMPI *);
    void initPhases(Params&, SmileiMPI *);
    void initParticles(Params&, std::vector<Species*>& vecSpecies);
    void initTestParticles(Params&, std::vector<Species*>&);
    
    //! field dump output
    unsigned int fieldDump_every;
    
    //! time-averaged field dump output
    unsigned int avgfieldDump_every;
    
    //! number of time-steps for time-averaging of fields
    unsigned int ntime_step_avg;
    
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
	
};

#endif
