// -----------------------------------------------------------------------------
//
//! \file ParticleFactory.h
//
//! \brief Factory to select the right Particles type
//
// -----------------------------------------------------------------------------

#ifndef PARTICLES_PARTICLESFACTORY_H
#define PARTICLES_PARTICLESFACTORY_H

#include "Params.h"
#include "Particles.h"

class ParticlesFactory
{
public:
    static Particles* create( const Params& parameters,
                              const Patch&  a_parent_patch );
};

#endif
