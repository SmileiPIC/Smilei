// -----------------------------------------------------------------------------
//
//! \file ParticleInjector.h
//
//! \brief Class with functions to inject particles
//
// -----------------------------------------------------------------------------

#ifndef PARTICLEINJECTOR_H
#define PARTICLEINJECTOR_H

#include <vector>
#include <string>

#include "Params.h"
#include "Profile.h"

//! class ParticleInjector
class ParticleInjector
{
public:
    
    // -----------------------------------------------------------------------------
    //  1. Constructor/destructor

    //! ParticleInjector constructor
    ParticleInjector( Params &, Patch * );

    //! ParticleInjector destructor
    virtual ~ParticleInjector();

    
    // -----------------------------------------------------------------------------
    //  2. Injector parameters
    
    //! number of the injector
    unsigned int injector_number_;

    //! Name of this injector
    std::string name_;

    //! kind/name of the species associated to this injector
    std::string species_name_;

    //! number of the species associated to this injector
    unsigned int species_number_;

    //! box side to where inject particles
    std::string box_side_;
    
    //! position initialization of the injector
    std::string position_initialization_;

    //! momentum initialization type, possible values: "cold" or "maxwell-juettner"
    std::string momentum_initialization_;

    //! Boolean to know if we initialize particles with positions of another injector
    bool position_initialization_on_injector_;
    
    //! Index of the species where position initialization is made
    int position_initialization_on_injector_index_;

    //! vector of velocity profiles (vx, vy, vz)
    std::vector<Profile *> velocity_profile_;
    
    //! vector of temperature profiles (Tx, Ty, Tz)
    std::vector<Profile *> temperature_profile_;

    //! Type of density profile ("nb" or "charge")
    std::string density_profile_type_;

    //! density profile
    Profile *density_profile_;

    //! time profile
    Profile *time_profile_;
    
    //! Profile for the particles per cell
    Profile *particles_per_cell_profile_;

    //! Pointer toward regular number of particles array
    std::vector<int> regular_number_array_;

    // -----------------------------------------------------------------------------
    //  3. Methods

    //! Method returning the number of the associated species
    inline unsigned int getSpeciesNumber() const
    {
        return species_number_;
    }

    //! Return if the injector is from Xmin
    inline bool isXmin()
    {
        return (box_side_ == "xmin");
    }

    //! Return if the injector is from Xmax
    inline bool isXmax()
    {
        return (box_side_ == "xmax");
    }

    //! Return if the injector is from Ymin
    inline bool isYmin()
    {
        return (box_side_ == "ymin");
    }

    //! Return if the injector is from Ymax
    inline bool isYmax()
    {
        return (box_side_ == "ymax");
    }

    //! Return if the injector is from Zmin
    inline bool isZmin()
    {
        return (box_side_ == "zmin");
    }

    //! Return if the injector is from Zmax
    inline bool isZmax()
    {
        return (box_side_ == "zmax");
    }

};

#endif
