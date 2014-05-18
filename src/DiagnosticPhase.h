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
    //! position of the particle
	std::vector<double> pos;
    //! momentum of the particle
	std::vector<double> mom;
    //! weight of the particle
	double weight;
    //! charge of the particle
	short charge;
};

//! this class holds what all kind of phase space diag should do
class DiagnosticPhase {

public:

    //! creator, each subclass will pick from phaseStructure what it needs
    DiagnosticPhase(phaseStructure);
    ~DiagnosticPhase(){};
	
    //! all diags should have this every parameter
	unsigned int every;
	
	//! this is to remember on which species calculate the phase space
	std::vector<std::string> my_species;
	
	//! this will update internal Field with the particle
	virtual void run(partStruct& my_part)=0;
	
	//! this will white the diagnostic header to the hdf5 file
	virtual void writeAttributes(hid_t gid)=0;
	
    //! this will write the internal Field to the file
	virtual void writeData(unsigned int timestep, hid_t gid)=0;

};
#endif
