#ifndef DiagnosticPhase_H
#define DiagnosticPhase_H

#include <hdf5.h>
#include "Tools.h"
#include "DiagParams.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

//! this strucs holds the basics of a particle
//!\todo check if this is slow (TV to JR)
struct partStruct {
	std::vector<double> pos;
	std::vector<double> mom;
	double weight;
	short charge;
};

//! this class holds what all kind of phase space diag should do
class DiagnosticPhase {

public:

    DiagnosticPhase(phaseStructure);
    ~DiagnosticPhase(){};
	
    //! all diags should have this every parameter
	unsigned int every;
	
	//! this is to remember on which species calculate the phase space
	std::vector<std::string> my_species;
	
	//! this will update internal Field with the particle
	virtual void doSomething(partStruct& my_part)=0;
	
	//! this will write the internal Field to the file
	virtual void writeData(unsigned int timestep, hid_t gid)=0;

};
#endif
