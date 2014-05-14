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
#include "DiagnosticProbe0D.h"
#include "DiagnosticProbe1D.h"
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
    Diagnostic(PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi, Interpolator* interp);
    //! destructor
    ~Diagnostic();
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
 
private:
    int num_CPUs;
    DiagnosticScalar diagScal;
    unsigned int everyScalar;

    DiagnosticProbe0D probe0D;

    DiagnosticProbe1D probe1D;
    
    Interpolator* interp_;

	DiagnosticPhaseSpace diagPhaseSpace;
	
};

#endif
