// -----------------------------------------------------------------------------
//
//! \file createParticles.h
//
//! \brief Class with functions to create particles
//
// -----------------------------------------------------------------------------

#ifndef CREATEPARTICLE_H
#define CREATEPARTICLE_H

#include <cstring>
#include <string>
#include <iostream>
#include <vector>

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "Field3D.h"
#include "ParticleCreatorContainer.h"

using namespace std;

class CreateParticles
{
public:
    
    //! Constructor
    CreateParticles() {};
    
    //! Destructor
    ~CreateParticles() {};
    
    //! Creation of the particle properties in the given particle vector `particles`
    static int create( struct particles_creator particle_creator,
                       Particles * particles,
                       Species * species,
                       vector<unsigned int> n_space_to_create,
                       Params &params,
                       Patch *patch,
                       int new_cell_idx );
    
    //! Creation of the particle positions
    static void createPosition( string position_initialization,
                              Particles * particles,
                              Species * species,
                              unsigned int nPart,
                              unsigned int iPart, double *indexes, Params &params );
    
    //! Creation of the particle momentum
    static void createMomentum( string momentum_initialization,
                            Particles * particles,
                            Species * species,
                            unsigned int nPart,
                            unsigned int iPart,
                            double *temp,
                            double *vel);
    
    //! Creation of the particle weight
    static void createWeight( Particles * particles,
                                               unsigned int nPart,
                                               unsigned int iPart,
                                               double n_real_particles );
    
    // For all particles in a mesh initialize its charge state
    static void createCharge( Particles * particles, Species * species,
                                       unsigned int nPart, unsigned int iPart, double q );
    
private:

    //! Provides a Maxwell-Juttner distribution of energies
    static vector<double> maxwellJuttner( Species * species, unsigned int npoints, double temperature );
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvF[1000];
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvH[1000];
};

#endif
