
#ifndef COLLISIONALFUSIONDD_H
#define COLLISIONALFUSIONDD_H

#include <vector>

#include "CollisionalNuclearReaction.h"

class CollisionalFusionDD : public CollisionalNuclearReaction
{

public:
    //! Constructor
    CollisionalFusionDD( Params&, std::vector<Species*>&, double );
    //! Cloning Constructor
    CollisionalFusionDD( CollisionalNuclearReaction * );
    //! Destructor
    ~CollisionalFusionDD() {};
    
    std::string name() override {
        return "Nuclear reaction: D-D fusion";
    };
    
    //! Method to apply the nuclear reaction
    double crossSection( double log_ekin ) override;
    //! Method to prepare the products of the reaction
    void makeProducts(  Random* random, double ekin, double log_ekin, double tot_charge, NuclearReactionProducts &products ) override;
    
    //! static parameters to read the database
    static const double a1, a2, a3, npointsm1;
    static const int npoints;
    static const double DB_log_crossSection[50];
    
private:
    unsigned int index_He_;
    unsigned int index_n_;
    
};

#endif
