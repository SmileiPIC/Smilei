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
    //! creator (just keep a copy of extfield_struct in my_struct)
    ExtFieldProfile(ExtFieldStructure &extfield_struct) : my_struct(extfield_struct) {};
    //! destructor (empty)
    virtual ~ExtFieldProfile() {};
    //! real operator for the profile (takes a position and returns the value of the imposed field)
    virtual double operator() (std::vector<double>)=0;
    
    
protected:
    //! copy of the structure that keeps the profile parameters
    ExtFieldStructure my_struct;
    
};//END class

#endif
