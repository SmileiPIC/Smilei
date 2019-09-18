#include "CollisionalNuclearReaction.h"

#include "Collisions.h"
#include "Species.h"
#include "Patch.h"

#include <cmath>


using namespace std;

// Constructor
CollisionalNuclearReaction::CollisionalNuclearReaction( Params *params, int product_species, Particles* particles )
{
    product_species_ = product_species;
    rate_multiplier_ = 1.;
    n_reactions_ = 0;
    if( params ) {
        product_particles.initialize( 0, *particles );
    }
}

// Cloning Constructor
CollisionalNuclearReaction::CollisionalNuclearReaction( CollisionalNuclearReaction *CNR )
{
    product_species_ = CNR->product_species_;
    rate_multiplier_ = CNR->rate_multiplier_;
    n_reactions_ = CNR->n_reactions_;
    product_particles.initialize( 0, CNR->product_particles );
}

// Finish the reaction
void CollisionalNuclearReaction::finish(
    Params &params, Patch *patch, std::vector<Diagnostic *> &localDiags,
    bool intra_collisions, vector<unsigned int> sg1, vector<unsigned int> sg2,
    double npairs, int itime
) {
    // Move new particles in place
    patch->vecSpecies[product_species_]->importParticles( params, patch, product_particles, localDiags );
    
    // Remove reactants that have been (very rare)
    for( unsigned int is = 0; is < sg1.size(); is++ ) {
        patch->vecSpecies[sg1[is]]->eraseWeightlessParticles();
    }
    if( ! intra_collisions ) {
        for( unsigned int is = 0; is < sg2.size(); is++ ) {
            patch->vecSpecies[sg2[is]]->eraseWeightlessParticles();
        }
    }
    
    // Update the rate multiplier
    double goal = (double) params.n_time * (double) n_reactions_ / ( (double) itime * npairs );
    //double goal = 1000000. * (double) n_reactions_ / ( (double) itime * npairs );
    if( goal > 2. ) {
        rate_multiplier_ *= 0.01;
    } else if( goal > 1.5 ) {
        rate_multiplier_ *= 0.1;
    } else if( goal > 1. ) {
        rate_multiplier_ *= 0.3;
    } else if( rate_multiplier_ < 1.e14 ) {
        if( goal < 0.01 ) {
            rate_multiplier_ *= 10.;
        } else if( goal < 0.1 ) {
            rate_multiplier_ *= 2.;
        } 
    }
}
