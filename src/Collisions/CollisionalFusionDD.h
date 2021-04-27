
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
    double crossSection( double log_ekin );
    //! Method to prepare the products of the reaction
    void makeProducts( double U, double etot, double log_ekin, double q, Particles *&p3, Particles *&p4, double &p3_COM, double &p4_COM, double &q3, double &q4, double &sinX, double &cosX ) override;
    
    std::string name() override { return "D-D fusion"; };
    
    //! static parameters to read the database
    static const double a1, a2, a3, npointsm1;
    static const int npoints;
    static const double DB_log_crossSection[50];
    
private:
    
    
};

#endif
