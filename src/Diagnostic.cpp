#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "DiagParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic( PicParams* params,  DiagParams* diagparams, SmileiMPI* smpi) :
scalars(params, diagparams, smpi),
probes(params, diagparams, smpi),
diagPhaseSpace(params, diagparams, smpi)
{
}

Diagnostic::~Diagnostic () {
    scalars.close();
    probes.close();
	diagPhaseSpace.close();
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp) {
    scalars.run(timestep, EMfields, vecSpecies);
    probes.run(timestep, EMfields, interp);
	diagPhaseSpace.run(timestep, vecSpecies);
}

