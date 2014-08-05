/*
 * Diagnostic.h
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#ifndef Diagnostic_H
#define Diagnostic_H

#include <vector>
#include <fstream>

#include "DiagnosticScalar.h"
#include "Interpolator.h"
#include "DiagnosticProbe.h"
#include "DiagnosticPhaseSpace.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;
class Species;


//! class holder for all the diagnostics: scalars, probes(0D, 1D, 2D and 3D) and phase-space
class Diagnostic {

public:
    //! creator called from main
    Diagnostic(PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi);
    
    //! destructor
    ~Diagnostic();
    
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp);
 
    //! get a particular scalar
    double getScalar(std::string name);
        
private:
    
    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace phases;
	
};

#endif
