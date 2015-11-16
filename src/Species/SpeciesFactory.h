#ifndef SPECIESFACTORY_H
#define SPECIESFACTORY_H

#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"

#include "Params.h"
#include "Patch.h"

#include "Tools.h"


class SpeciesFactory {
public:
    static Species* create(Params& params, int ispec, Patch* patch) {
        Species* sp = NULL;
        if (params.species_param[ispec].dynamics_type=="norm") {
            // Species with Boris dynamics
	    sp = new Species_norm(params, ispec, patch);

        } else if (params.species_param[ispec].dynamics_type=="rrll") {
            // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
            sp = new Species_rrll(params, ispec, patch);
        } // END if

	if (params.species_param[ispec].isTest) {
	    int locNbrParticles = sp->getNbrOfParticles();
	    int sz(1);
	    MPI_Comm_size( MPI_COMM_WORLD, &sz );
	    int* allNbrParticles = new int[sz];
	    MESSAGE( "Internal synchro to do to gather per by patch master" );
	    MPI_Gather( &locNbrParticles, 1, MPI_INTEGER, allNbrParticles, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
	    int nParticles(0);

	    if ( patch->isMaster() ) {
		nParticles =  allNbrParticles[0];
		for (int irk=1 ; irk<sz ; irk++){
		    allNbrParticles[irk] += nParticles;
		    nParticles = allNbrParticles[irk];
		}
		for (int irk=sz-1 ; irk>0 ; irk--){
		    allNbrParticles[irk] = allNbrParticles[irk-1];
		}
		allNbrParticles[0] = 0;

	    }

	    int offset(0);
	    MPI_Scatter(allNbrParticles, 1 , MPI_INTEGER, &offset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD );
	    sp->particles->addIdOffsets(offset);
	}

        return sp;
    }

    static std::vector<Species*> createVector(Params& params, Patch* patch) {
        std::vector<Species*> vecSpecies;
        vecSpecies.resize(params.species_param.size());

        Species *electron_species=NULL;

        // create species
        unsigned int nPart;
        for (unsigned int ispec=0 ; ispec<params.species_param.size() ; ispec++) {
	    vecSpecies[ispec] = SpeciesFactory::create(params, ispec, patch);
            if (params.species_param[ispec].species_type=="electron") {
                electron_species=vecSpecies[ispec];
            }
            nPart = vecSpecies[ispec]->getNbrOfParticles();
        } // END for ispec

        // add the found electron species to the ionizable species
        for (unsigned int ispec=0 ; ispec<params.species_param.size() ; ispec++) {
            if (vecSpecies[ispec]->Ionize)  {
                if (electron_species) {
                    vecSpecies[ispec]->electron_species=electron_species;
                    //Compilation fails with "pointer to incomplete class type is not allowed" if following line is uncommented.
                    //PMESSAGE(2,patch->Hindex(),"Added electron species to species " << vecSpecies[ispec]->species_param.species_type);
                } else {
                    ERROR("Ionization needs a species called \"electron\" to be defined");
                }
            }
        } // END for ispec

        return vecSpecies;
    }

};

#endif
