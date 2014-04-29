#ifndef DiagnosticPhase2D_H
#define DiagnosticPhase2D_H

#include "DiagnosticPhase.h"
#include "Field2D.h"

class PicParams;
class DiagParams;
class ElectroMagn;

//! this class holds all the phase projections that can be represented as 2d matrix
class DiagnosticPhase2D : public DiagnosticPhase {

public:
    DiagnosticPhase2D(phaseStructure phaseStruct);

	// specialized method to write a field 2d
	void writeData(unsigned int timestep, hid_t gid);
	
	//! by now this is the easiest way 2d classes holds field2d and they know how to write it
	Field2D my_data;
	
protected:
	//! first component of the phasespace min and max
	double firstmin,firstmax;
	//! second component of the phasespace min and max
	double secondmin,secondmax;
	//! number of bins for the first and second component of the 2d phasespace
	unsigned int firstnum,secondnum;
	
};
#endif
