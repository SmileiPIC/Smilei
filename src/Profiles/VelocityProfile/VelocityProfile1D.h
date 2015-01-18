#ifndef VelocityProfile1D_H
#define VelocityProfile1D_H


#include "VelocityProfile.h"
#include <cmath>

//  --------------------------------------------------------------------------------------------------------------------
//! 1D density profile class
//  --------------------------------------------------------------------------------------------------------------------
class VelocityProfile1D : public VelocityProfile
{
    
public:
    VelocityProfile1D(ProfileSpecies &params);
    ~VelocityProfile1D() {};
    double operator() (std::vector<double>);
    
private:
    
    
};//END class

#endif
