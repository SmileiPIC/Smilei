#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"
#include "DiagnosticParticles.h"

using namespace std;

Diagnostic::Diagnostic( PicParams &picParams, DiagParams &dParams , SmileiMPI* smpi) :
scalars(picParams, dParams, smpi),
probes(picParams, dParams, smpi),
phases(picParams, dParams, smpi)
{
    dtimer[0].init(smpi, "scalars");
    dtimer[1].init(smpi, "probes");
    dtimer[2].init(smpi, "phases");

}

Diagnostic::~Diagnostic () {

    MESSAGE(0, "Time in scalars : " << dtimer[0].getTime() );
    MESSAGE(0, "Time in probes : " << dtimer[1].getTime() );
    MESSAGE(0, "Time in phases : " << dtimer[2].getTime() );

    scalars.close();
    probes.close();
	phases.close();
	DiagnosticParticles::closeAll();
}

double Diagnostic::getScalar(string name){
    return scalars.getScalar(name);
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp, SmileiMPI *smpi) {
    dtimer[0].restart();
    scalars.run(timestep, EMfields, vecSpecies, smpi);
    dtimer[0].update();

    dtimer[1].restart();
    probes.run(timestep, EMfields, interp);
    dtimer[1].update();

    dtimer[2].restart();
    phases.run(timestep, vecSpecies);
    dtimer[2].update();

    DiagnosticParticles::runAll(timestep, vecSpecies, smpi);

}

