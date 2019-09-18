
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
    CollisionalNuclearReaction( Params *, int, Particles* );
    //! Cloning Constructor
    CollisionalNuclearReaction( CollisionalNuclearReaction * );
    //! Destructor
    virtual ~CollisionalNuclearReaction() {};
    
    //! Method to prepare the nuclear reaction
    virtual void prepare() {};
    //! Method to test the occurence of the nuclear reaction
    virtual bool occurs( double U, double coeff, double m1, double m2, double g1, double g2, double &etot, double &log_ekin, double &W ) = 0;
    //! Method to prepare the products of the reaction
    virtual void makeProducts( double U, double ekin, double log_ekin, double q, Particles *&p3, Particles *&p4, double &p3_COM, double &p4_COM, double &q3, double &q4, double &cosX ) = 0;
    //! Method to finish the nuclear reaction and put new electrons in place
    virtual void finish( Params &, Patch *, std::vector<Diagnostic *> &, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime );
    
    //! New electrons temporary species
    Particles product_particles;
    
    //! Coefficient to adapt the rate of create of new particles. Automatically adjusted
    double rate_multiplier_;
    
    //! Number of succesfull reactions since the beginning (in this patch)
    unsigned int n_reactions_;
    
private:
    
    //! Species where products are sent
    int product_species_;
    
};

//! Class to make empty reaction objects
class CollisionalNoNuclearReaction : public CollisionalNuclearReaction
{
public:
    CollisionalNoNuclearReaction() : CollisionalNuclearReaction( NULL, -1, NULL ) {};
    ~CollisionalNoNuclearReaction() {};
    
    void prepare() override {};
    bool occurs( double U, double coeff, double m1, double m2, double g1, double g2, double &etot, double &log_ekin, double &W ) override {
        return false;
    };
    void makeProducts( double U, double ekin, double log_ekin, double q, Particles *&p3, Particles *&p4, double &p3_COM, double &p4_COM, double &q3, double &q4, double &cosX ) override {};
    void finish( Params &params, Patch *patch, std::vector<Diagnostic *> &diags, bool, std::vector<unsigned int>, std::vector<unsigned int>, double npairs, int itime ) override {};
};


#endif
