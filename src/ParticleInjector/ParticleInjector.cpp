#include "ParticleInjector.h"

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for ParticleInjector
// ---------------------------------------------------------------------------------------------------------------------
ParticleInjector::ParticleInjector( Params &params, Patch *patch ):
    velocity_profile_( 3, NULL ),
    temperature_profile_( 3, NULL )
{
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destructor for ParticleInjector
// ---------------------------------------------------------------------------------------------------------------------
ParticleInjector::~ParticleInjector() {

    for( unsigned int i=0; i<velocity_profile_.size(); i++ ) {
        delete velocity_profile_[i];
    }
    for( unsigned int i=0; i<temperature_profile_.size(); i++ ) {
        delete temperature_profile_[i];
    }
}
