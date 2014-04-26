#ifndef DiagnosticPhase_H
#define DiagnosticPhase_H

#include <hdf5.h>
#include "Tools.h"
#include "DiagParams.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

struct partStruct {
	std::vector<double> pos;
	std::vector<double> mom;
	double weight;
	short charge;
};

class DiagnosticPhase {

public:

    DiagnosticPhase(phaseStructure);
    ~DiagnosticPhase();
	
	unsigned int every;
	std::vector<std::string> my_species;
	
	//! this will update internal Field with the particle
	virtual void doSomething(partStruct& my_part)=0;
	
	//! this will write the internal Field to the file
	virtual void writeData(unsigned int timestep, SmileiMPI* smpi, hid_t gid)=0;

};
#endif
