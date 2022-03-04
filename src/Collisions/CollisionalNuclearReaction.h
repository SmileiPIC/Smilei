#ifndef COLLISIONALNUCLEARREACTION_H
#define COLLISIONALNUCLEARREACTION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"
#include "Random.h"
#include "NuclearReactionProducts.h"

class Patch;

class CollisionalNuclearReaction
{

public:
    //! Constructor
    CollisionalNuclearReaction( Params *, std::vector<Species *>*, double );
    //! Cloning Constructor
    CollisionalNuclearReaction( CollisionalNuclearReaction * );
    //! Destructor
    virtual ~CollisionalNuclearReaction();
    
    //! Prepare the nuclear reaction
    virtual void prepare() {
        tot_probability_ = 0.;
    };
    //! Calculates the cross-section vs energy
    virtual double crossSection( double log_ekin ) = 0;
    //! Test the occurence of the nuclear reaction
    virtual bool occurs( double U, double coeff, double etot, double &log_ekin, double &W );
    //! Prepare the products of the reaction
    virtual void makeProducts( Random* random, double ekin, double log_ekin, double tot_charge, NuclearReactionProducts &products ) = 0;
    //! Finish the nuclear reaction and put new electrons in place
    virtual void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime );
    
    //! Get the name of the reaction
    virtual std::string name() = 0;
    
    //! Temporary species for reaction products
    std::vector<Particles *> product_particles_;
    
    //! Coefficient to adapt the rate of create of new particles
    double rate_multiplier_;
    
    //! True if rate multiplier isn't automatically adjusted
    double auto_multiplier_;
    
    //! sum of probabilities in this patch
    double tot_probability_;
    
private:
    
    //! Species where products are sent
    std::vector<unsigned int> product_ispecies_;
    
    //! Mass of reaction products
    std::vector<double> product_mass_;
    
};

//! Class to make empty reaction objects
class CollisionalNoNuclearReaction : public CollisionalNuclearReaction
{
public:
    CollisionalNoNuclearReaction() : CollisionalNuclearReaction( NULL, NULL, 0. ) {};
    ~CollisionalNoNuclearReaction() {};
    
    void prepare() override {};
    double crossSection( double log_ekin ) override { return 0.; };
    bool occurs( double U, double coeff, double etot, double &log_ekin, double &W ) override {
        return false;
    };
    void makeProducts( Random* random, double ekin, double log_ekin, double tot_charge, NuclearReactionProducts &products ) override {};
    void finish( Params &params, Patch *patch, std::vector<Diagnostic *> &diags, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime ) override {};
    std::string name() override { return ""; }
};


#endif
