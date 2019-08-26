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

// ---------------------------------------------------------------------------------------------------------------------
//! Return the particle_creator structure from the properties of the injector
// ---------------------------------------------------------------------------------------------------------------------
void ParticleInjector::getParticlesCreator() {
    
    particles_creator_.position_initialization_ = position_initialization_;
    particles_creator_.position_initialization_on_species_ = position_initialization_on_injector_;
    particles_creator_.momentum_initialization_ = momentum_initialization_;
    particles_creator_.velocity_profile_.resize(velocity_profile_.size());
    for (unsigned int i = 0 ; i < velocity_profile_.size() ; i++) {
        particles_creator_.velocity_profile_[i] = new Profile(velocity_profile_[i]);
    }
    particles_creator_.temperature_profile_.resize(temperature_profile_.size());
    for (unsigned int i = 0 ; i < temperature_profile_.size() ; i++) {
        particles_creator_.temperature_profile_[i] = new Profile(temperature_profile_[i]);
    }
    particles_creator_.density_profile_ = new Profile(density_profile_);
    particles_creator_.density_profile_type_ = density_profile_type_;
    particles_creator_.time_profile_ = new Profile(time_profile_);
    particles_creator_.particles_per_cell_profile_ = new Profile(particles_per_cell_profile_);

}
