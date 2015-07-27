#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic( PicParams &picParams, DiagParams &dParams , SmileiMPI* smpi, Patch* patch) :
scalars(picParams, dParams, patch),
probes(picParams, dParams, patch),
phases(picParams, dParams, smpi)
{
}

Diagnostic::~Diagnostic () {
    scalars.close();
    probes.close();
	phases.close();
}

double Diagnostic::getScalar(string name){
    return scalars.getScalar(name);
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp) {
    scalars.run(timestep, EMfields, vecSpecies);
    probes.run(timestep, EMfields, interp);
	phases.run(timestep, vecSpecies);
}

