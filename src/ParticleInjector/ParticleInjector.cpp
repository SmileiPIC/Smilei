#include "ParticleInjector.h"

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for ParticleInjector
// ---------------------------------------------------------------------------------------------------------------------
ParticleInjector::ParticleInjector( Params &params, Patch *patch ):
    velocity_profile_( 3, NULL ),
    temperature_profile_( 3, NULL ),
    density_profile_( NULL ),
    time_profile_( NULL ),
    particles_per_cell_profile_( NULL )
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
    if( particles_per_cell_profile_ ) {
        delete particles_per_cell_profile_;
    }
    if( time_profile_ ) {
        delete time_profile_;
    }
    if( density_profile_ ) {
        delete density_profile_;
    }
}
