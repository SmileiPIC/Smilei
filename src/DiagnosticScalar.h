#ifndef DIAGNOSTICSCALAR_H
#define DIAGNOSTICSCALAR_H

#include "Tools.h"
//#include "Diagnostic.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <hdf5.h>

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;
class Species;



class DiagnosticScalar{

public:
	
	DiagnosticScalar(SmileiMPI* smpi);
	//destructor
	~DiagnosticScalar(){};	
	
	void run(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	void compute(int timestep, ElectroMagn* EMfields, std::vector<Species*>&);
	void write(int timestep,std::vector<Species*>&);
		
private:
	std::ofstream fout;
	SmileiMPI* smpi_;
	int num_CPUs;
	MPI_Datatype mpi_data_scalar;
	std::vector<double> vecScalar;
	
};

#endif
