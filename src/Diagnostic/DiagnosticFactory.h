#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "DiagnosticParticles.h"
#include "DiagnosticProbes.h"
#include "DiagnosticScalar.h"
#include "DiagnosticTrack.h"


class DiagnosticFactory {
public:
    /*static Diagnostic* create(Params& params, SmileiMPI* smpi, unsigned int  ipatch) {
	return 
    }*/

    static std::vector<Diagnostic*> createGlobalDiagnostics(Params& params, SmileiMPI* smpi, Patch* patch) {

	std::vector<Diagnostic*> vecDiagnostics;

	vecDiagnostics.push_back( new DiagnosticScalar(params, smpi, patch, 0) ); // 0 : 1 scalar only (useless)

	for (unsigned int n_diag_particles = 0; n_diag_particles < PyTools::nComponents("DiagParticles"); n_diag_particles++) {
	    // append new diagnostic object to the list
	    vecDiagnostics.push_back( new DiagnosticParticles(params, smpi, patch, n_diag_particles) );
	}
        
	return vecDiagnostics;

    } // END createGlobalDiagnostics


static std::vector<Diagnostic*> createLocalDiagnostics(Params& params, SmileiMPI* smpi, Patch* patch) {

	std::vector<Diagnostic*> vecDiagnostics;

	for (unsigned int n_diag_probes = 0; n_diag_probes < PyTools::nComponents("DiagProbe"); n_diag_probes++) {
	    vecDiagnostics.push_back( new DiagnosticProbes(params, smpi, patch, n_diag_probes) );
	}

	// writable particles initialization
	// loop species and make a new diag if particles have to be dumped
	for(unsigned int trackIdx=0; trackIdx<patch->vecSpecies.size(); trackIdx++) {
	    if ( patch->vecSpecies[trackIdx]->particles->tracked ) {
	      vecDiagnostics.push_back( new DiagnosticTrack(params, smpi, patch, vecDiagnostics.size() ) ); // trackIdx not used, no python parsing to init
	    }
	}

	return vecDiagnostics;

    } // END createLocalDiagnostics

};

#endif

