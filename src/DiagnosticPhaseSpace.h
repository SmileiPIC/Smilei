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
	SmileiMPI *smpi_;
	std::vector<DiagnosticPhase*> vecDiagPhase;
	
	hid_t fileId;

	unsigned int ndim;
};
#endif
