#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic( PicParams &picParams, InputData &ifile , SmileiMPI* smpi) :
params(picParams,ifile),
scalars(picParams, params, smpi),
probes(picParams, params, smpi),
phases(picParams, params, smpi)
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

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp, SmileiMPI *smpi) {
    scalars.run(timestep, EMfields, vecSpecies, smpi);
    probes.run(timestep, EMfields, interp);
	phases.run(timestep, vecSpecies);
}

