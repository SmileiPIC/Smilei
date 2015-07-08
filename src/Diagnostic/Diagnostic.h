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
#include "Timer.h"

class PicParams;
class SmileiMPI;
class ElectroMagn;
class Species;


//! class holder for all the diagnostics: scalars, probes(0D, 1D, 2D and 3D) and phase-space
class Diagnostic {

public:
    //! creator called from main
    Diagnostic(PicParams&, InputData&, SmileiMPI *smpi);
    
    //! destructor
    ~Diagnostic(){};
    
    //! print timers
    void printTimers(SmileiMPI *smpi, double tottime);
    
    //! close all diags
    void closeAll();

    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp, SmileiMPI *smpi);
 
    //! get a particular scalar
    double getScalar(std::string name);

    std::vector<Timer> dtimer;
    
    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace phases;
	
    std::vector<DiagnosticParticles*> vecDiagnosticParticles;
    
    DiagParams params;
};

#endif
