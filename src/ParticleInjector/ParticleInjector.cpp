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

// ---------------------------------------------------------------------------------------------------------------------
//! Return the particle_creator structure from the properties of the injector
// ---------------------------------------------------------------------------------------------------------------------
struct particles_creator ParticleInjector::getParticlesCreator() {

    struct particles_creator particles_creator;
    
    particles_creator.position_initialization_ = position_initialization_;
    particles_creator.position_initialization_on_species_ = position_initialization_on_injector_;
    particles_creator.momentum_initialization_ = momentum_initialization_;
    particles_creator.velocity_profile_.resize(velocity_profile_.size());
    for (unsigned int i = 0 ; i < velocity_profile_.size() ; i++) {
        particles_creator.velocity_profile_[i] = new Profile(velocity_profile_[i]);
    }
    particles_creator.temperature_profile_.resize(temperature_profile_.size());
    for (unsigned int i = 0 ; i < temperature_profile_.size() ; i++) {
        particles_creator.temperature_profile_[i] = new Profile(temperature_profile_[i]);
    }
    particles_creator.density_profile_ = new Profile(density_profile_);
    particles_creator.density_profile_type_ = density_profile_type_;
    particles_creator.time_profile_ = new Profile(time_profile_);
    
    return particles_creator;
}
