#include "ElectroMagnBC.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "LaserParams.h"
#include "LaserProfile.h"
#include "Tools.h"

using namespace std;

// Constructor for ElectromagnBC
ElectroMagnBC::ElectroMagnBC( PicParams &params, LaserParams &laser_params )
{
    // check for laser conditions
    laser_.resize(laser_params.n_laser);

    for (unsigned int i=0; i<laser_.size(); i++) {
        DEBUG(5,"Initializing Laser "<<i);        
        laser_[i] = new LaserProfile(params,laser_params, i);
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

