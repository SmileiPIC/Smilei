#ifndef DIAGFACTORY_H
#define DIAGFACTORY_H

#include "DiagParticles.h"
#include "DiagProbes.h"
#include "DiagScalar.h"
#include "DiagTrack.h"


class DiagFactory {
public:
    /*static Diag* create(Params& params, SmileiMPI* smpi, unsigned int  ipatch) {
	return 
    }*/

    static std::vector<Diag*> createGlobalDiags(Params& params, SmileiMPI* smpi, Patch* patch) {

	std::vector<Diag*> vecDiags;

	vecDiags.push_back( new DiagScalar(params, smpi, patch, 0) ); // 0 : 1 scalar only (useless)

	for (unsigned int n_diag_particles = 0; n_diag_particles < PyTools::nComponents("DiagParticles"); n_diag_particles++) {
	    // append new diagnostic object to the list
	    vecDiags.push_back( new DiagParticles(params, smpi, patch, n_diag_particles) );
	}
        
	return vecDiags;

    } // END createGlobalDiags


static std::vector<Diag*> createLocalDiags(Params& params, SmileiMPI* smpi, Patch* patch) {

	std::vector<Diag*> vecDiags;

	for (unsigned int n_diag_probes = 0; n_diag_probes < PyTools::nComponents("DiagProbe"); n_diag_probes++) {
	    vecDiags.push_back( new DiagProbes(params, smpi, patch, n_diag_probes) );
	}

	// writable particles initialization
	// loop species and make a new diag if particles have to be dumped
	for(unsigned int trackIdx=0; trackIdx<patch->vecSpecies.size(); trackIdx++) {
	    if ( patch->vecSpecies[trackIdx]->particles->tracked ) {
		vecDiags.push_back( new DiagTrack(params, smpi, patch, trackIdx ) ); // trackIdx not used, no python parsing to init
	    }
	}

	return vecDiags;

    } // END createLocalDiags

};

#endif

