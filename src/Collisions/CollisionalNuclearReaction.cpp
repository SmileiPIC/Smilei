#include "CollisionalNuclearReaction.h"

#include "Collisions.h"
#include "Species.h"
#include "Particles.h"
#include "Patch.h"

#include <cmath>


using namespace std;

// Constructor
CollisionalNuclearReaction::CollisionalNuclearReaction(
    Params &params,
    vector<Species *> &product_species,
    double rate_multiplier
    )
{
    if( rate_multiplier == 0. ) {
        auto_multiplier_ = true;
        rate_multiplier_ = 1.;
    } else {
        auto_multiplier_ = false;
        rate_multiplier_ = rate_multiplier;
    }
    product_particles_.resize(0);
    product_ispecies_.resize(0);
    product_mass_.resize(0);
    for( unsigned int i=0; i<product_species.size(); i++ ) {
        if( product_species[i] != NULL ) {
            product_particles_.push_back( new Particles() );
            product_particles_.back()->initialize( 0, *( product_species[i]->particles ) );
            product_ispecies_.push_back( product_species[i]->species_number_ );
            product_mass_.push_back( product_species[i]->mass_ );
        }
    }
    coeff2_ = 2.817940327e-15*params.reference_angular_frequency_SI/299792458.; // re omega / c
}

// Cloning Constructor
CollisionalNuclearReaction::CollisionalNuclearReaction( CollisionalNuclearReaction *CNR )
{
    product_ispecies_ = CNR->product_ispecies_;
    rate_multiplier_ = CNR->rate_multiplier_;
    auto_multiplier_ = CNR->auto_multiplier_;
    product_particles_.resize( CNR->product_particles_.size(), NULL );
    for( unsigned int i=0; i<CNR->product_particles_.size(); i++ ) {
        product_particles_[i] = new Particles();
        product_particles_[i]->initialize( 0, *CNR->product_particles_[i] );
    }
    coeff2_ = CNR->coeff2_;
}

// Cloning Constructor
CollisionalNuclearReaction::~CollisionalNuclearReaction()
{
    for( unsigned int i=0; i<product_particles_.size(); i++ ) {
        delete product_particles_[i];
    }
    product_particles_.clear();
}

void CollisionalNuclearReaction::apply( Random *random, BinaryProcessData &D )
{
    double ekin = D.m1 * (D.gamma1_COM-1.) + D.m2 * (D.gamma2_COM-1.);
    double log_ekin = log( ekin );
    
    // Interpolate the total cross-section at some value of ekin = m1(g1-1) + m2(g2-1)
    double cs = crossSection( log_ekin );
    
    // Calculate probability for reaction
    double prob = coeff2_ * D.vrel_corr * D.dt_correction * cs * rate_multiplier_;
    tot_probability_ += prob;
    npairs_tot_ ++;
    if( random->uniform() > exp( -prob ) ) {
        
        // Reaction occurs
        
        double W = D.minW / rate_multiplier_;
        
        // Reduce the weight of both reactants
        // If becomes zero, then the particle will be discarded later
        D.p1->weight( D.i1 ) -= W;
        D.p2->weight( D.i2 ) -= W;
        D.minW -= 0.;
        D.maxW -= 0.;
        
        // Get the magnitude and the angle of the outgoing products in the COM frame
        NuclearReactionProducts products;
        double tot_charge = D.p1->charge( D.i1 ) + D.p2->charge( D.i2 );
        makeProducts( random, ekin, log_ekin, tot_charge, products );
        
        // Calculate new weights
        double newW1, newW2;
        if( tot_charge != 0. ) {
            double weight_factor = W / tot_charge;
            newW1 = D.p1->charge( D.i1 ) * weight_factor;
            newW2 = D.p2->charge( D.i2 ) * weight_factor;
        } else {
            newW1 = W;
            newW2 = 0.;
        }
        
        // For each product
        double p_perp = sqrt( D.px_COM*D.px_COM + D.py_COM*D.py_COM );
        double newpx_COM=0, newpy_COM=0, newpz_COM=0;
        for( unsigned int iproduct=0; iproduct<products.particles.size(); iproduct++ ){
            // Calculate the deflection in the COM frame
            if( iproduct < products.cosPhi.size() ) { // do not recalculate if all products have same axis
                if( p_perp > 1.e-10*D.p_COM ) { // make sure p_perp is not too small
                    double inv_p_perp = 1./p_perp;
                    newpx_COM = ( D.px_COM * D.pz_COM * products.cosPhi[iproduct] - D.py_COM * D.p_COM * products.sinPhi[iproduct] ) * inv_p_perp;
                    newpy_COM = ( D.py_COM * D.pz_COM * products.cosPhi[iproduct] + D.px_COM * D.p_COM * products.sinPhi[iproduct] ) * inv_p_perp;
                    newpz_COM = -p_perp * products.cosPhi[iproduct];
                } else { // if p_perp is too small, we use the limit px->0, py=0
                    newpx_COM = D.p_COM * products.cosPhi[iproduct];
                    newpy_COM = D.p_COM * products.sinPhi[iproduct];
                    newpz_COM = 0.;
                }
                // Calculate the deflection in the COM frame
                newpx_COM = newpx_COM * products.sinX[iproduct] + D.px_COM *products.cosX[iproduct];
                newpy_COM = newpy_COM * products.sinX[iproduct] + D.py_COM *products.cosX[iproduct];
                newpz_COM = newpz_COM * products.sinX[iproduct] + D.pz_COM *products.cosX[iproduct];
            }
            // Go back to the lab frame and store the results in the particle array
            double vcp = D.COM_vx * newpx_COM + D.COM_vy * newpy_COM + D.COM_vz * newpz_COM;
            double momentum_ratio = products.new_p_COM[iproduct] / D.p_COM;
            double term6 = momentum_ratio*D.term1*vcp + sqrt( products.new_p_COM[iproduct]*products.new_p_COM[iproduct] + 1. ) * D.COM_gamma;
            double newpx = momentum_ratio * newpx_COM + D.COM_vx * term6;
            double newpy = momentum_ratio * newpy_COM + D.COM_vy * term6;
            double newpz = momentum_ratio * newpz_COM + D.COM_vz * term6;
            // Make new particle at position of particle 1
            if( newW1 > 0. ) {
                products.particles[iproduct]->makeParticleAt( *D.p1, D.i1, newW1, products.q[iproduct], newpx, newpy, newpz );
            }
            // Make new particle at position of particle 2
            if( newW2 > 0. ) {
                products.particles[iproduct]->makeParticleAt( *D.p2, D.i2, newW2, products.q[iproduct], newpx, newpy, newpz );
            }
        }
        
    } // end nuclear reaction

}


// Finish the reaction
void CollisionalNuclearReaction::finish(
    Params &params, Patch *patch, std::vector<Diagnostic *> &localDiags,
    bool intra_collisions, vector<unsigned int> sg1, vector<unsigned int> sg2, int itime
) {
    // Move new particles in place
    for( unsigned int i=0; i<product_particles_.size(); i++ ) {
        patch->vecSpecies[product_ispecies_[i]]->importParticles( params, patch, *product_particles_[i], localDiags );
    }
    
    // Remove reactants that have fully reacted (very rare)
    for( unsigned int is = 0; is < sg1.size(); is++ ) {
        patch->vecSpecies[sg1[is]]->eraseWeightlessParticles();
    }
    if( ! intra_collisions ) {
        for( unsigned int is = 0; is < sg2.size(); is++ ) {
            patch->vecSpecies[sg2[is]]->eraseWeightlessParticles();
        }
    }
    
    if( auto_multiplier_ ) {
        // Update the rate multiplier
        double goal = tot_probability_ * (double) params.n_time / npairs_tot_;
        if( goal > 0. ) {
            if( goal > 2. ) {
                rate_multiplier_ *= 0.01;
            } else if( goal > 1.5 ) {
                rate_multiplier_ *= 0.1;
            } else if( goal > 1. ) {
                rate_multiplier_ *= 0.3;
            } else if( rate_multiplier_ < 1.e14 ) {
                if( goal < 0.01 ) {
                    rate_multiplier_ *= 8.;
                } else if( goal < 0.1 ) {
                    rate_multiplier_ *= 2.;
                }
            }
            if( rate_multiplier_ < 1. ) {
                rate_multiplier_ = 1.;
            }
        }
    }
}
