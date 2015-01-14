#ifndef ExtFieldProfile1D_H
#define ExtFieldProfile1D_H


#include "ExtFieldProfile.h"

//  --------------------------------------------------------------------------------------------------------------------
//! 1D ExtField profile class
//  --------------------------------------------------------------------------------------------------------------------
class ExtFieldProfile1D : public ExtFieldProfile
{
    
public:
    ExtFieldProfile1D(ExtFieldStructure &extfield_struct);
    ~ExtFieldProfile1D() {};
    double operator() (std::vector<double>);
    
private:
    
    
};//END class

#endif
