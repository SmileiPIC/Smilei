
#ifndef COLLISIONALFUSIONDD_H
#define COLLISIONALFUSIONDD_H

#include <vector>

#include "CollisionalNuclearReaction.h"

class CollisionalFusionDD : public CollisionalNuclearReaction
{

public:
    //! Constructor
    CollisionalFusionDD( Params*, std::vector<Particles*>, std::vector<unsigned int>, double );
    //! Cloning Constructor
    CollisionalFusionDD( CollisionalNuclearReaction * );
    //! Destructor
    ~CollisionalFusionDD() {};
    
    //! Method to apply the nuclear reaction
    double crossSection( double log_ekin ) override;
    //! Method to prepare the products of the reaction
    void makeProducts( Random* random, double ekin, double log_ekin, double tot_charge, std::vector<Particles *> &particles, std::vector<double> &p_COM, std::vector<short> &q, std::vector<double> &sinX, std::vector<double> &cosX ) override;
    
    std::string name() override { return "D-D fusion"; };
    
    //! static parameters to read the database
    static const double a1, a2, a3, npointsm1;
    static const int npoints;
    static const double DB_log_crossSection[50];
    
private:
    
    
};

#endif
