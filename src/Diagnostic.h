/*
 * Diagnostic.h
 *
 *  Created on: 3 juil. 2013
 *      Author: jderouil
 */

#ifndef Diagnostic_H
#define Diagnostic_H

#include <vector>
#include <iostream>
#include <fstream>
#include <hdf5.h>

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;
class Species;

class Diagnostic {
public:
	Diagnostic(PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi);
	~Diagnostic(){};
	void compute(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	
	
	void write(int timestep);

private:
	SmileiMPI* smpi_;
	std::vector<double> data_;
	int every;
	std::ofstream fout;
	std:: vector<double> mean_values;
	std::vector<unsigned int> mean_weight;
	
	int num_CPUs;
	
};

#endif
