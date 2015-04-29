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
    Diagnostic(SmileiMPI* smpi);
    
    //! destructor
    ~Diagnostic(){};
    
    //! close all diags
    void closeAll();
    
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp, SmileiMPI *smpi);
 
    //! get a particular scalar
    double getScalar(std::string name);

    Timer dtimer[4];

    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace phases;
	
    std::vector<DiagnosticParticles*> vecDiagnosticParticles;
    
};

#endif
