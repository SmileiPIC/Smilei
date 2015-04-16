#ifndef VelocityProfile_H
#define VelocityProfile_H

#include <cmath>
#include "PicParams.h"



//  --------------------------------------------------------------------------------------------------------------------
//! Class VelocityProfile
//  --------------------------------------------------------------------------------------------------------------------
class VelocityProfile
{
public:
    VelocityProfile(ProfileSpecies &my_params) : prof_params(my_params) {
        if (prof_params.profile.empty()) {
           // prof_params.profile="constant";
            prof_params.profile="none";
        }        
    };
    virtual ~VelocityProfile() {};
    virtual double operator() (std::vector<double>)=0;
        
protected:
    ProfileSpecies prof_params;
    
};//END class

#endif
