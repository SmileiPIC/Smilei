#ifndef COLLISIONALNUCLEARREACTION_H
#define COLLISIONALNUCLEARREACTION_H

#include <vector>

#include "Tools.h"
#include "Species.h"
#include "Params.h"

class Patch;

class CollisionalNuclearReaction
{

public:
    //! Constructor
    CollisionalNuclearReaction( Params *, std::vector<Particles*>*, std::vector<unsigned int>*, double );
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
    virtual bool occurs( double U, double coeff, double m1, double m2, double g1, double g2, double &etot, double &log_ekin, double &W );
    //! Prepare the products of the reaction
    virtual void makeProducts( double U, double ekin, double log_ekin, double q, Particles *&p3, Particles *&p4, double &p3_COM, double &p4_COM, double &q3, double &q4, double &cosX ) = 0;
    //! Finish the nuclear reaction and put new electrons in place
    virtual void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime );
    
    //! Get the name of the reaction
    virtual std::string name() = 0;
    
    //! New electrons temporary species
    std::vector<Particles *> product_particles_;
    
    //! Coefficient to adapt the rate of create of new particles
    double rate_multiplier_;
    
    //! True if rate multiplier isn't automatically adjusted
    double auto_multiplier_;
    
    //! sum of probabilities in this patch
    double tot_probability_;
    
private:
    
    //! Species where products are sent
    std::vector<unsigned int> product_species_;
    
};

//! Class to make empty reaction objects
class CollisionalNoNuclearReaction : public CollisionalNuclearReaction
{
public:
    CollisionalNoNuclearReaction() : CollisionalNuclearReaction( NULL, NULL, NULL, 0. ) {};
    ~CollisionalNoNuclearReaction() {};
    
    void prepare() override {};
    double crossSection( double log_ekin ) override { return 0.; };
    bool occurs( double U, double coeff, double m1, double m2, double g1, double g2, double &etot, double &log_ekin, double &W ) override {
        return false;
    };
    void makeProducts( double U, double ekin, double log_ekin, double q, Particles *&p3, Particles *&p4, double &p3_COM, double &p4_COM, double &q3, double &q4, double &cosX ) override {};
    void finish( Params &params, Patch *patch, std::vector<Diagnostic *> &diags, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime ) override {};
    std::string name() override { return ""; }
};


#endif
