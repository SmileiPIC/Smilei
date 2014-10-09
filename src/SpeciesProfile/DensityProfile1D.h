#ifndef DensityProfile1D_H
#define DensityProfile1D_H


#include "DensityProfile.h"

//  --------------------------------------------------------------------------------------------------------------------
//! 1D density profile class
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile1D : public DensityProfile
{
    
public:
    DensityProfile1D(SpeciesStructure &params);
    ~DensityProfile1D() {};
    double operator() (std::vector<double>);
    
private:
    
    
};//END class

#endif
