#ifndef DensityProfile1D_H
#define DensityProfile1D_H


#include "DensityProfile.h"

//  --------------------------------------------------------------------------------------------------------------------
//! lkjskljaslkdjaskl
//  --------------------------------------------------------------------------------------------------------------------
class DensityProfile1D : public DensityProfile
{
    
public:
    DensityProfile1D(){};
    ~DensityProfile1D() {};
    double operator() (PicParams*, unsigned int, std::vector<double>);
    
private:
    
    
};//END class

#endif
