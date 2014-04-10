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
#include "DiagnosticProbe0D.h"
#include "DiagnosticPhase1D.h"

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
    //! destructor (empty)
    ~Diagnostic() {};
    //! check if at timestep diagnostics must be called
    void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
    void closeAll();


private:
    int num_CPUs;
    DiagnosticScalar diagScal;
    unsigned int everyScalar;

    DiagnosticProbe0D probe0D;
    unsigned int everyProbe0D;
    Interpolator* interp_;

	DiagnosticPhase1D phase1D;
	unsigned int everyPhase1D;
	
    unsigned int everyMap;
};

#endif
