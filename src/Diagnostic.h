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
	void runAllDiags (int timestep, ElectroMagn* EMfields, std::vector<Species*>&);


private:
	unsigned int everyScalar;
	int num_CPUs;
	DiagnosticScalar DiagScal;
	
};

#endif
