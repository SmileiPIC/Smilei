#ifndef DiagnosticPhase_H
#define DiagnosticPhase_H

#include <hdf5.h>
#include "Tools.h"
#include "DiagParams.h"

class PicParams;
class SmileiMPI;
class DiagParams;
class ElectroMagn;

class DiagnosticPhase {

public:

    DiagnosticPhase(phaseStructure, hid_t gid);
    ~DiagnosticPhase();
	
	void close();
	unsigned int every;
	hid_t groupID;
	std::vector<std::string> my_species;
};
#endif
