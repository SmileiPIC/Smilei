#ifndef TemperatureProfile_H
#define TemperatureProfile_H

#include <cmath>
#include "PicParams.h"



//  --------------------------------------------------------------------------------------------------------------------
//! Class TemperatureProfile
//  --------------------------------------------------------------------------------------------------------------------
class TemperatureProfile
{
public:
    TemperatureProfile(ProfileSpecies &my_params) : prof_params(my_params) {
        if (prof_params.profile.empty()) {
            //prof_params.profile="constant";
            prof_params.profile="none";
        }        
    };
    virtual ~TemperatureProfile() {};
    virtual double operator() (std::vector<double>)=0;
        
protected:
    ProfileSpecies prof_params;
    
};//END class

#endif
