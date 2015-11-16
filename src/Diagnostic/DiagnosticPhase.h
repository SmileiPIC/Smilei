#ifndef DiagnosticPhase_H
#define DiagnosticPhase_H

#include <hdf5.h>
#include "Tools.h"
#include "Field2D.h"

class Params;
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
    
    inline double lor_fact() { 
        return sqrt(1.0+pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2));
    }
};

struct phaseStructure {
	//!string defining the kind oh phase projections
    std::vector<std::string> kind;
    
    //! phase output every (every phase diagnostic must have this)
    unsigned int every;
    
    //! phase output from tmin
    double tmin;
    
    //! phase output to tmin
    double tmax;
    
    //! compression level using zlib [0-9] (0 deactvate compression)
    unsigned int deflate;
    
	//! vector of pointer to species on which the phase diag will be applied (if omitted, it will be for all)
	std::vector<std::string> species;
	
    //! minimum position
	std::vector<double> pos_min;
    //! max position
	std::vector<double> pos_max;
    //! number of positions
	std::vector <unsigned int> pos_num;
    
    //! minimum momentum
	std::vector<double> mom_min;
    //! max momentum
	std::vector<double> mom_max;
    //! number of momenta
	std::vector <unsigned int> mom_num;
	
    //! minimum Lorentz factor
	std::vector<double> lor_min;
    //! max Lorentz factor
	std::vector<double> lor_max;
    //! number of Lorentz factors
	std::vector <unsigned int> lor_num;
	
};


//! this class holds all the phase projections that can be represented as 2d matrix
class DiagnosticPhase {

public:
    //! creator
    DiagnosticPhase(phaseStructure phaseStruct);
    ~DiagnosticPhase();

    //! this will write the internal Field to the file
	void writeData();
	
    //! all diags should have this every parameter
	unsigned int every;
	
	//! this is to remember on which species calculate the phase space
	std::vector<std::string> my_species;
	
	//! this will update internal Field with the particle
	virtual void run(partStruct& my_part)=0;
	
	//! by now this is the easiest way 2d classes holds field2d and they know how to write it
	Field2D my_data;
	
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
	
    //! HDF identifier of data
    hid_t dataId;

    //! dt not passed to this, tmin/tmax do not run for this for now
    double tmin;
    double tmax;
    double dt;

};
#endif
