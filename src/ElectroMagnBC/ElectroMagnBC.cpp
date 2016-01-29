
#include <cstdlib>
#include <iostream>
#include <string>

#include "ElectroMagnBC.h"

#include "Params.h"
#include "LaserParams.h"
#include "LaserProfile.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;

// Constructor for ElectromagnBC
ElectroMagnBC::ElectroMagnBC( Params &params, LaserParams &laser_params, Patch* patch )
{
    // check for laser conditions
    laser_.resize(laser_params.laser_param.size());
    
    for (unsigned int i=0; i<laser_.size(); i++) {
        DEBUG("Initializing Laser "<<i);        
        laser_[i] = new LaserProfile(params,laser_params, i, patch);
    }

    // time step
    dt = params.timestep;
}

// Destructor for ElectromagnBC
ElectroMagnBC::~ElectroMagnBC()
{
    for (unsigned int i=0; i< laser_.size(); i++) {
        delete laser_[i];
    }
}

// Disable all lasers when using moving window
void ElectroMagnBC::laserDisabled()
{
    for (unsigned int i=0; i< laser_.size(); i++) {
        laser_[i]->disabled();
    }
}

