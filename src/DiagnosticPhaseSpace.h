#ifndef DiagnosticPhaseSpace_H
#define DiagnosticPhaseSpace_H

#include <cmath>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <hdf5.h>

#include "Tools.h"

#include "Species.h"
#include "Interpolator.h"
#include "Particles.h"

#include "DiagnosticPhase.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

class DiagnosticPhaseSpace {

public:

    DiagnosticPhaseSpace(PicParams* params, DiagParams* diagParams, SmileiMPI* smpi);
    ~DiagnosticPhaseSpace();

	void run(int timestep, std::vector<Species*>& vecSpecies);
	
	void close();
private:
    
    //! this vector will hold all the diagnostics created
	std::vector<DiagnosticPhase*> vecDiagPhase;
	
    //! this is the hdf5 file id (we need to keep it to close at the right time)
	hid_t fileId;

    //! this is always handy to know (number of particle dimension)
	unsigned int ndim;
	
    //! this holds in which hdf5 groupID (hid_t) will the data for each species (string) of DiagnosticPhase is written
    std::map<DiagnosticPhase*, std::map<std::string,hid_t> >mapGroupId; 
};
#endif
