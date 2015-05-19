#ifndef TemperatureProfile1D_H
#define TemperatureProfile1D_H


#include "TemperatureProfile.h"
#include <cmath>

//  --------------------------------------------------------------------------------------------------------------------
//! 1D density profile class
//  --------------------------------------------------------------------------------------------------------------------
class TemperatureProfile1D : public TemperatureProfile
{
    
public:
    TemperatureProfile1D(ProfileSpecies &params);
    ~TemperatureProfile1D() {};
    double operator() (std::vector<double>);
    
private:
    
    
};//END class

#endif
