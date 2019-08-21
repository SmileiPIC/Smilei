// -----------------------------------------------------------------------------
//
//! \file ParticleCreatorContainer.h
//
//! \brief Class with a container structure to pass the particle creator properties
//
// -----------------------------------------------------------------------------

#ifndef PARTICLECREATORCONTAINER_H
#define PARTICLECREATORCONTAINER_H

#include <cstring>
#include <vector>
#include "Profile.h"

//! Structure that contains particle properties for their creation
struct particles_creator {
    // Initialization with the positions of another species
    bool position_initialization_on_species_;
    // Position initialization type
    std::string position_initialization_;
    // Momentum initialization type
    std::string momentum_initialization_;
    //! vector of velocity profiles (vx, vy, vz)
    std::vector<Profile *> velocity_profile_;
    //! vector of temperature profiles (Tx, Ty, Tz)
    std::vector<Profile *> temperature_profile_;
    //! Density profile
    Profile * density_profile_;
    //! Type of profile
    std::string density_profile_type_;
};

#endif
