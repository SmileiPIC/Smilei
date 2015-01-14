#ifndef ExtFieldProfile_H
#define ExtFieldProfile_H

#include <cmath>
#include "ExtFieldParams.h"



//  --------------------------------------------------------------------------------------------------------------------
//! Class ExtFieldProfile
//  --------------------------------------------------------------------------------------------------------------------
class ExtFieldProfile
{
public:
    ExtFieldProfile(ExtFieldStructure &extfield_struct) : my_struct(extfield_struct) {};
    virtual ~ExtFieldProfile() {};
    virtual double operator() (std::vector<double>)=0;
    
    
protected:
    ExtFieldStructure my_struct;
    
};//END class

#endif
