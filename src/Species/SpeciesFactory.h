#ifndef SPECIESFACTORY_H
#define SPECIESFACTORY_H

#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Patch.h"

#include "Tools.h"

class SpeciesFactory {
public:
  static Species* create(PicParams& params, int ispec, SmileiMPI* smpi, Patch* patch) {
        Species* sp = NULL;
        if (params.species_param[ispec].dynamics_type=="norm") {
            // Species with Boris dynamics
	    sp = new Species_norm(params, ispec, smpi, patch);

        } else if (params.species_param[ispec].dynamics_type=="rrll") {
            // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
            sp = new Species_rrll(params, ispec, smpi, patch);
        } // END if

        return sp;
    }

  static std::vector<Species*> createVector(PicParams& params, SmileiMPI* smpi, Patch* patch) {
        std::vector<Species*> vecSpecies;
        vecSpecies.resize(params.n_species);

        Species *electron_species=NULL;

        // create species
        unsigned int nPart;
        for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
	  vecSpecies[ispec] = SpeciesFactory::create(params, ispec, smpi, patch);
            if (params.species_param[ispec].species_type=="electron") {
                electron_species=vecSpecies[ispec];
            }
            nPart = vecSpecies[ispec]->getNbrOfParticles();
        } // END for ispec

        // add the found electron species to the ionizable species
        for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
            if (vecSpecies[ispec]->Ionize)  {
                if (electron_species) {
                    vecSpecies[ispec]->electron_species=electron_species;
                    PMESSAGE(2,smpi->getRank(),"Added electron species to species " << vecSpecies[ispec]->species_param.species_type);
                } else {
                    ERROR("Ionization needs a species called \"electron\" to be defined");
                }
            }
        } // END for ispec

        return vecSpecies;
    }

};

#endif
