#ifndef COLLISIONALNUCLEARREACTION_H
#define COLLISIONALNUCLEARREACTION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"
#include "Random.h"
#include "BinaryProcess.h"
#include "NuclearReactionProducts.h"

class Patch;

class CollisionalNuclearReaction : public BinaryProcess
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
    void apply( Random *random, BinaryProcessData &D );
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
