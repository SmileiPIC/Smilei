#ifndef DensityProfile_H
#define DensityProfile_H

#include <cmath>
#include "PicParams.h"
#include "pyHelper.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class DensityProfile
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile
{
public:
    DensityProfile(SpeciesStructure &params) : species_param(params) {};
    virtual ~DensityProfile() {};
    virtual double operator() (std::vector<double>)=0;
    
    
protected:
    SpeciesStructure species_param;
    
};//END class

#endif
