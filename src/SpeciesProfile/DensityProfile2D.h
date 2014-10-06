#ifndef DensityProfile2D_H
#define DensityProfile2D_H


#include "DensityProfile.h"

//  --------------------------------------------------------------------------------------------------------------------
//! lkjskljaslkdjaskl
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile2D : public DensityProfile
{
    
public:
    DensityProfile2D(SpeciesStructure &params) : DensityProfile(params) {};
    ~DensityProfile2D() {};
    double operator() (std::vector<double>);
    
    
private:
    
    
};//END class

#endif
