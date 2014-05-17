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


//! class holder of the two type of diagnostics scalars and map
class Diagnostic {

public:
    //! creator called from main
    Diagnostic(PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi);
    //! destructor
    ~Diagnostic();
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator *interp);
 
private:

    DiagnosticScalar scalars;

    DiagnosticProbe probes;

	DiagnosticPhaseSpace diagPhaseSpace;
	
};

#endif
