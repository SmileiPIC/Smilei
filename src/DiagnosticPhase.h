#ifndef DiagnosticPhase_H
#define DiagnosticPhase_H

#include <hdf5.h>
#include "Tools.h"
#include "DiagParams.h"
#include "Field2D.h"

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

//! this class holds all the phase projections that can be represented as 2d matrix
class DiagnosticPhase {

public:
    //! creator
    DiagnosticPhase(phaseStructure phaseStruct);
    ~DiagnosticPhase(){};

    //! this will write the internal Field to the file
	void writeData(unsigned int timestep, hid_t gid);
	
    //! all diags should have this every parameter
	unsigned int every;
	
	//! this is to remember on which species calculate the phase space
	std::vector<std::string> my_species;
	
	//! this will update internal Field with the particle
	virtual void run(partStruct& my_part)=0;
	
	//! this will white the diagnostic header to the hdf5 file
    void writeAttributes(hid_t did);

	//! by now this is the easiest way 2d classes holds field2d and they know how to write it
	Field2D my_data;
	
protected:
	//! first component of the phasespace min
	double firstmin;
	//! first component of the phasespace max
    double firstmax;
	//! number of bins for the first component of the 2d phasespace
	unsigned int firstnum;
	//! second component of the phasespace min
	double secondmin;
	//! second component of the phasespace max
    double secondmax;
	//! number of bins for the second component of the 2d phasespace
	unsigned int secondnum;
	
};
#endif
