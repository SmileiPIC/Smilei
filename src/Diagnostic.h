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
#include "DiagParams.h"

class PicParams;
class SmileiMPI;
class ElectroMagn;
class Species;


//! class holder for all the diagnostics: scalars, probes(0D, 1D, 2D and 3D) and phase-space
class Diagnostic {

public:
    //! creator called from main
    Diagnostic(PicParams &params,  InputData &ifile, SmileiMPI* smpi);
    
    //! destructor
    ~Diagnostic();
    
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp, SmileiMPI *smpi);
 
    //! get a particular scalar
    double getScalar(std::string name);
        
    DiagParams params;
private:
    
    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace phases;
	
};

#endif
