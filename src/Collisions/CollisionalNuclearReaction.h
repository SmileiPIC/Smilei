#ifndef COLLISIONALNUCLEARREACTION_H
#define COLLISIONALNUCLEARREACTION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"
#include "Random.h"
#include "NuclearReactionProducts.h"
#include "BinaryProcessData.h"

class Patch;

class CollisionalNuclearReaction
{

public:
    //! Constructor
    CollisionalNuclearReaction( Params &, std::vector<Species *>&, double );
    //! Cloning Constructor
    CollisionalNuclearReaction( CollisionalNuclearReaction * );
    //! Destructor
    virtual ~CollisionalNuclearReaction();
    
    //! Prepare the nuclear reaction
    void prepare() {
        tot_probability_ = 0.;
        npairs_tot_ = 0;
    };

    SMILEI_ACCELERATOR_DECLARE_ROUTINE
    void apply( Random *random, BinaryProcessData &D, uint32_t n )
    {
    // Not ported to GPU yet
    #ifndef SMILEI_ACCELERATOR_GPU
        for( uint32_t i = 0; i<n; i++ ) {
            
            double ekin = D.m[0][i] * D.gamma_tot_COM[i] - D.m[0][i] - D.m[1][i];
            double log_ekin = log( ekin );
            
            // Interpolate the total cross-section at some value of ekin = m1(g1-1) + m2(g2-1)
            double cs = crossSection( log_ekin );
            
            // Calculate probability for reaction
            double vrel_corr = D.p_COM[i] * D.gamma_tot_COM[i] / ( D.R[i]* D.gamma[0][i] * D.gamma[1][i] );
            double prob = coeff2_ * vrel_corr * D.dt_correction[i] * cs * rate_multiplier_;
            tot_probability_ += prob;
            
            if( random->uniform() > exp( -prob ) ) {
                
                // Reaction occurs
                
                double W = std::min( D.W[0][i], D.W[1][i] ) / rate_multiplier_;
                
                // Reduce the weight of both reactants
                // If becomes zero, then the particle will be discarded later
                D.W[0][i] -= W;
                D.W[1][i] -= W;
                
                // Get the magnitude and the angle of the outgoing products in the COM frame
                NuclearReactionProducts products;
                double tot_charge = D.q[0][i] + D.q[1][i];
                makeProducts( random, ekin, log_ekin, tot_charge, products );
                
                // Calculate new weights
                double newW1, newW2;
                if( tot_charge != 0. ) {
                    double weight_factor = W / tot_charge;
                    newW1 = D.q[0][i] * weight_factor;
                    newW2 = D.q[1][i] * weight_factor;
                } else {
                    newW1 = W;
                    newW2 = 0.;
                }
                
                // For each product
                double p_perp = sqrt( D.px_COM[i]*D.px_COM[i] + D.py_COM[i]*D.py_COM[i] );
                double newpx_COM=0, newpy_COM=0, newpz_COM=0;
                for( unsigned int iproduct=0; iproduct<products.particles.size(); iproduct++ ){
                    // Calculate the deflection in the COM frame
                    if( iproduct < products.cosPhi.size() ) { // do not recalculate if all products have same axis
                        if( p_perp > 1.e-10*D.p_COM[i] ) { // make sure p_perp is not too small
                            double inv_p_perp = 1./p_perp;
                            newpx_COM = ( D.px_COM[i] * D.pz_COM[i] * products.cosPhi[iproduct] - D.py_COM[i] * D.p_COM[i] * products.sinPhi[iproduct] ) * inv_p_perp;
                            newpy_COM = ( D.py_COM[i] * D.pz_COM[i] * products.cosPhi[iproduct] + D.px_COM[i] * D.p_COM[i] * products.sinPhi[iproduct] ) * inv_p_perp;
                            newpz_COM = -p_perp * products.cosPhi[iproduct];
                        } else { // if p_perp is too small, we use the limit px->0, py=0
                            newpx_COM = D.p_COM[i] * products.cosPhi[iproduct];
                            newpy_COM = D.p_COM[i] * products.sinPhi[iproduct];
                            newpz_COM = 0.;
                        }
                        // Calculate the deflection in the COM frame
                        newpx_COM = newpx_COM * products.sinX[iproduct] + D.px_COM[i] *products.cosX[iproduct];
                        newpy_COM = newpy_COM * products.sinX[iproduct] + D.py_COM[i] *products.cosX[iproduct];
                        newpz_COM = newpz_COM * products.sinX[iproduct] + D.pz_COM[i] *products.cosX[iproduct];
                    }
                    // Rescale by the new momentum magnitude
                    double momentum_ratio = products.new_p_COM[iproduct] / D.p_COM[i];
                    double new_gamma_COM = sqrt( products.new_p_COM[iproduct]*products.new_p_COM[iproduct] + 1. );
                    newpx_COM *= momentum_ratio;
                    newpy_COM *= momentum_ratio;
                    newpz_COM *= momentum_ratio;
                    // Go back to the lab frame and store the results in the particle array
                    double pp = ( D.px_tot[i] * newpx_COM + D.py_tot[i] * newpy_COM + D.pz_tot[i] * newpz_COM ) / ( D.gamma_tot[i] + D.gamma_tot_COM[i] );
                    double f = ( new_gamma_COM + pp ) / D.gamma_tot_COM[i];
                    double newpx = newpx_COM + f * D.px_tot[i];
                    double newpy = newpy_COM + f * D.py_tot[i];
                    double newpz = newpz_COM + f * D.pz_tot[i];
                    // Make new particle at position of particle 1
                    if( newW1 > 0. ) {
                        products.particles[iproduct]->makeParticleAt( *D.p[0][i], D.i[0][i], newW1, products.q[iproduct], newpx, newpy, newpz );
                    }
                    // Make new particle at position of particle 2
                    if( newW2 > 0. ) {
                        products.particles[iproduct]->makeParticleAt( *D.p[1][i], D.i[1][i], newW2, products.q[iproduct], newpx, newpy, newpz );
                    }
                }
                
            } // end nuclear reaction
        }
        
        npairs_tot_ += n;
    #endif
    }
    SMILEI_ACCELERATOR_DECLARE_ROUTINE_END
    
    void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool intra, std::vector<unsigned int> sg1, std::vector<unsigned int> sg2, int itime );
    virtual std::string name() = 0;
    
    //! Calculates the cross-section vs energy
    virtual double crossSection( double log_ekin ) = 0;
    //! Prepare the products of the reaction
    virtual void makeProducts( Random* random, double ekin, double log_ekin, double tot_charge, NuclearReactionProducts &products ) = 0;
    
    //! Temporary species for reaction products
    std::vector<Particles *> product_particles_;
    
    //! Coefficient to adapt the rate of create of new particles
    double rate_multiplier_;
    
    //! True if rate multiplier isn't automatically adjusted
    double auto_multiplier_;
    
    //! sum of probabilities in this patch
    double tot_probability_;
    
    unsigned int npairs_tot_;
private:
    
    //! Species where products are sent
    std::vector<unsigned int> product_ispecies_;
    
    //! Mass of reaction products
    std::vector<double> product_mass_;
    
    double coeff2_;
    
};


#endif
