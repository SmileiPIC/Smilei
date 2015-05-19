#ifndef DensityProfile2D_H
#define DensityProfile2D_H


#include "DensityProfile.h"

//  --------------------------------------------------------------------------------------------------------------------
//! 2D density profile class
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile2D : public DensityProfile
{
    
public:
    DensityProfile2D(SpeciesStructure &params);
    ~DensityProfile2D() {};
    double operator() (std::vector<double>);
    
    
private:
    
    
};//END class

#endif
