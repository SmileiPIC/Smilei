#ifndef TemperatureProfile2D_H
#define TemperatureProfile2D_H


#include "TemperatureProfile.h"
#include <cmath>

//  --------------------------------------------------------------------------------------------------------------------
//! 2D density profile class
//  --------------------------------------------------------------------------------------------------------------------
class TemperatureProfile2D : public TemperatureProfile
{
    
public:
    TemperatureProfile2D(ProfileSpecies &params);
    ~TemperatureProfile2D() {};
    double operator() (std::vector<double>);
    
private:
    
    
};//END class

#endif
