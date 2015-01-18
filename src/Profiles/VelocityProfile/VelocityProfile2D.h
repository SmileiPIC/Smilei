#ifndef VelocityProfile2D_H
#define VelocityProfile2D_H


#include "VelocityProfile.h"
#include <cmath>

//  --------------------------------------------------------------------------------------------------------------------
//! lkjskljaslkdjaskl
//  --------------------------------------------------------------------------------------------------------------------
class VelocityProfile2D : public VelocityProfile
{
    
public:
    VelocityProfile2D(ProfileSpecies &params);
    ~VelocityProfile2D() {};
    double operator() (std::vector<double>);
    
    
private:
    
    
};//END class

#endif
