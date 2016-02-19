#include "Diagnostic.h"

#include <string>
#include <iomanip>

#include <hdf5.h>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic(Params& params, Patch* patch, SmileiMPI* smpi) :
dtimer(5),
scalars(params, patch),
probes(params, patch),
phases(params, patch)
{
    
    // defining default values & reading diagnostic every-parameter
    // ------------------------------------------------------------
    print_every=params.n_time/10;
    PyTools::extract("print_every", print_every);
    
    if (!PyTools::extract("fieldDump_every", fieldDump_every)) {
        fieldDump_every=params.global_every;
        MESSAGE(1,"Activating all fields to dump");
    }
    
    avgfieldDump_every=params.res_time*10;
    if (!PyTools::extract("avgfieldDump_every", avgfieldDump_every)) avgfieldDump_every=params.global_every;
    
    //!\todo Define default behaviour : 0 or params.res_time
    //ntime_step_avg=params.res_time;
    ntime_step_avg=0;
    PyTools::extract("ntime_step_avg", ntime_step_avg);
    
    // scalars initialization
    //dtimer[0].init(smpi, "scalars");
    
    // probes initialization
    //dtimer[1].init(smpi, "probes");
    
    // phasespaces timer initialization
    //dtimer[2].init(smpi, "phases");
    
    // particles initialization
    //dtimer[3].init(smpi, "particles");
    for (unsigned int n_diag_particles = 0; n_diag_particles < PyTools::nComponents("DiagParticles"); n_diag_particles++) {
        // append new diagnostic object to the list
        vecDiagnosticParticles.push_back(new DiagnosticParticles(n_diag_particles, params, patch, patch->vecSpecies));
    }
        
    // writable particles initialization
    //dtimer[4].init(smpi, "track particles");
    // loop species and make a new diag if particles have to be dumped
    for(unsigned int i=0; i<patch->vecSpecies.size(); i++) {
        if (patch->vecSpecies[i]->particles->track_every > 0) {
            vecDiagnosticTrackParticles.push_back(new DiagnosticTrackParticles(params, patch, patch->vecSpecies[i]));
        }
    }
}

void Diagnostic::closeAll (Patch* patch) {
    
    scalars.closeFile();
    probes.close();
    phases.close();
    
}

void Diagnostic::printTimers (Patch *patch, double tottime) {
    
    double coverage(0.);
    if ( patch->isMaster() ) {
        for (unsigned int i=0 ; i<dtimer.size() ; i++) {
            coverage += dtimer[i].getTime();
        }
    }
    MESSAGE(0, "\n Time in diagnostics : \t"<< tottime <<"\t" << coverage/tottime*100. << "% coverage" );
    if ( patch->isMaster() ) {
        for (unsigned int i=0 ; i<dtimer.size() ; i++) {
            dtimer[i].print(tottime) ;
        }
    }
}

double Diagnostic::getScalar(string name){
    return scalars.getScalar(name);
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp) {
    //dtimer[0].restart();
    scalars.run(timestep, EMfields, vecSpecies);
    //dtimer[0].update();
    
    //dtimer[1].restart();
    probes.run(timestep, EMfields, interp);
    //dtimer[1].update();
    
    //dtimer[2].restart();
    phases.run(timestep, vecSpecies);
    //dtimer[2].update();
    
    // run all the particle diagnostics
    //dtimer[3].restart();
    for (unsigned int i=0; i<vecDiagnosticParticles.size(); i++)
        vecDiagnosticParticles[i]->run(timestep, vecSpecies);
    //dtimer[3].update();

    // run all the track particle diagnostics
    //dtimer[4].restart();
    for (unsigned int i=0; i<vecDiagnosticTrackParticles.size(); i++)
        vecDiagnosticTrackParticles[i]->run(timestep);
    //dtimer[4].update();
    
}
