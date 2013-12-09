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
#include <hdf5.h>
#include "DiagnosticScalar.h"
<<<<<<< HEAD
=======
#include "Interpolator.h"
>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
#include "Probe0D.h"

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
	//! destructor (empty)
	~Diagnostic(){};
	//! check if at timestep diagnostics must be called
	void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&, Interpolator* interp);


private:
	int num_CPUs;
	DiagnosticScalar diagScal;
    unsigned int everyScalar;
<<<<<<< HEAD
=======
	
    Probe0D probe0d;
	unsigned int everyProbe0D;
    
    unsigned int everyMap;
>>>>>>> d9430c98fb22fc65390629334be5bfeebc514228
	
    Probe0D probe0d;
	unsigned int everyProbe0D;
    
    Interpolator* interp_;

    unsigned int everyMap;

};

#endif
