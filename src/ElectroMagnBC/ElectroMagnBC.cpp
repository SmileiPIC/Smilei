
#include <cstdlib>
#include <iostream>
#include <string>

#include "ElectroMagnBC.h"

#include "Params.h"
#include "Laser.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;

// Constructor for ElectromagnBC
ElectroMagnBC::ElectroMagnBC( Params &params, Patch *patch, unsigned int _min_max ) :
    min_max( _min_max )
{
    vecLaser.resize( 0 );
    
    // time step
    dt = params.timestep;
}

// Destructor for ElectromagnBC
ElectroMagnBC::~ElectroMagnBC()
{
    for( unsigned int i=0; i< vecLaser.size(); i++ ) {
        delete vecLaser[i];
    }
    vecLaser.clear();
}


// Disable all lasers when using moving window
void ElectroMagnBC::laserDisabled()
{
    //for (unsigned int i=0; i< vecLaser.size(); i++) {
    //    vecLaser[i]->disable();
    //}
    for( unsigned int i=0; i< vecLaser.size(); i++ ) {
        delete vecLaser[i];
    }
    vecLaser.resize( 0 );
}

