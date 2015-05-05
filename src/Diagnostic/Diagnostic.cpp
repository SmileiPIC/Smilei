#include "Diagnostic.h"

#include <string>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic(SmileiMPI* smpi) :
scalars(smpi),
probes(smpi),
phases(smpi)
{
    dtimer[0].init(smpi, "scalars");
    dtimer[1].init(smpi, "probes");
    dtimer[2].init(smpi, "phases");
    dtimer[3].init(smpi, "particles");

}

void Diagnostic::closeAll () {

    MESSAGE(0, "Time in diags : ");
    MESSAGE(0, "scalars : " << dtimer[0].getTime() );
    MESSAGE(0, "probes : " << dtimer[1].getTime() );
    MESSAGE(0, "phases : " << dtimer[2].getTime() );
    MESSAGE(0, "particles : " << dtimer[3].getTime() );

    scalars.close();
    probes.close();
	phases.close();
    
    for (int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->close();
    
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
    
    // run all the particle diagnostics
    dtimer[3].restart();
    for (int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->run(timestep, vecSpecies, smpi);
    dtimer[3].update();

}

