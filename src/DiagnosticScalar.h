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
	void compute(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	void init(ElectroMagn* EMfields, std::vector<Species*>&);
	void write(int timestep,std::vector<Species*>&);
		
private:
	std::ofstream fout;
	SmileiMPI* smpi_;
	std::vector<std::vector<spec_scalar_data> > mpi_spec_scalar_data;
};

#endif
