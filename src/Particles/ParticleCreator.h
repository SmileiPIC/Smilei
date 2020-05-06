// -----------------------------------------------------------------------------
//
//! \file ParticleCreator.h
//
//! \brief Class with functions to create particles
//
// -----------------------------------------------------------------------------

#ifndef PARTICLECREATOR_H
#define PARTICLECREATOR_H

#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "ParticleInjector.h"
#include "Field3D.h"

class ParticleCreator
{
public:
    
    //! Constructor
    ParticleCreator();
    
    //! Destructor
    ~ParticleCreator();
    
    //! Associate this particle creator object to the specified particle injector
    void associate( ParticleInjector * particle_injector, Particles * particles, Species * species);
    
    //! Associate this particle creator object to the specified species
    void associate( Species * species );
    
    //! Creation of the particle properties in the given particle vector `particles`
    int create( std::vector<unsigned int> n_space_to_create,
                       Params &params,
                       Patch *patch,
                       int new_cell_idx,
                       unsigned int itime );
    
    //! Creation of the particle positions
    static void createPosition( std::string position_initialization,
                              Particles * particles,
                              Species * species,
                              unsigned int nPart,
                              unsigned int iPart, double *indexes, Params &params );
    
    //! Creation of the particle momentum
    static void createMomentum( std::string momentum_initialization,
                            Particles * particles,
                            Species * species,
                            unsigned int nPart,
                            unsigned int iPart,
                            double *temp,
                            double *vel);
    
    //! Creation of the particle weight
    static void createWeight( std::string position_initialization,
                            Particles * particles,
                            unsigned int nPart,
                            unsigned int iPart,
                            double n_real_particles,
                            Params &params );

    //! Regulate of the particle weight with particle position
    static void regulateWeightwithPositionAM( Particles * particles, std::string position_initialization_on_species_type_, double dr );
    
    // For all particles in a mesh initialize its charge state
    static void createCharge( Particles * particles, Species * species,
                                       unsigned int nPart, unsigned int iPart, double q );
    
    //! Pointer the associated particle vector
    Particles * particles_;
    
    //! Pointer to the associated species
    Species * species_;
    
    //! Initialization with the positions of another species
    bool position_initialization_on_species_;

    //! Initialization type of the species where position initialization is made
    std::string position_initialization_on_species_type_;
    
    //! Flag if initialized in particles of a species
    bool initialized_in_species_;
    
    //! Position initialization type
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

    //! Time profile
    Profile * time_profile_;

    //! Particles per cell
    Profile * particles_per_cell_profile_;
    
    //! Flag for the addition of the energy coming from the created particles
    bool add_new_particle_energy_;
    
private:

    //! Provides a Maxwell-Juttner distribution of energies
    static std::vector<double> maxwellJuttner( Species * species, unsigned int npoints, double temperature );
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvF[1000];
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvH[1000];
};

#endif
