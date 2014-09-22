#ifndef DensityProfile_H
#define DensityProfile_H

#include <cmath>
#include "PicParams.h"



//  --------------------------------------------------------------------------------------------------------------------
//! Class DensityProfile
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile
{
public:
    DensityProfile() {};
    virtual ~DensityProfile() {};
    virtual double operator() (PicParams*, unsigned int, std::vector<double>)=0;
    
    
private:
    
};//END class

#endif
