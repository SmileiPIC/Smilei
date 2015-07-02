#include "Diagnostic.h"

#include <string>
#include <iomanip>

#include <hdf5.h>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic(PicParams& params, InputData &ifile, SmileiMPI *smpi) :
dtimer(4),
scalars(smpi),
probes(smpi),
phases(smpi),
params(this, params,ifile,smpi)
{
    dtimer[0].init(smpi, "scalars");
    dtimer[1].init(smpi, "probes");
    dtimer[2].init(smpi, "phases");
    dtimer[3].init(smpi, "particles");
    
}

void Diagnostic::closeAll () {
    
    scalars.close();
    probes.close();
    phases.close();
    
    for (int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->close();
    
}

void Diagnostic::printTimers (SmileiMPI *smpi, double tottime) {
    
    double coverage(0.);
    MESSAGE(0, "Time in diagnostics : \t"<< tottime);
    if ( smpi->isMaster() )
        for (int i=0 ; i<dtimer.size() ; i++) {
            dtimer[i].print(tottime) ;
            coverage += dtimer[i].getTime();
        }
    MESSAGE(0, "\t" << setw(12) << "Coverage\t\t" << coverage/tottime*100. << " %" );    
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

