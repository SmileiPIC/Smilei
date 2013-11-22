#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include "Tools.h"
//#include "Diagnostic.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>
#include "Species.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;



class DiagnosticScalar{

public:
	
	DiagnosticScalar(PicParams* params, SmileiMPI* smpi);
	//destructor
	~DiagnosticScalar(){};	
	
	void run(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	void compute_gather(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	void write(int timestep,std::vector<Species*>&);
		
private:
	std::ofstream fout;
	SmileiMPI* smpi_;
	//! scalars [iCpu][iSpec]
	std::vector<std::vector<std::map<std::string, double> > > mpi_spec_scalars;
};

#endif
