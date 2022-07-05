// -----------------------------------------------------------------------------
//
//! \file ParticleFactory.h
//
//! \brief Factory to select the right Particles type
//
// -----------------------------------------------------------------------------

#ifndef PARTICLESFACTORY_H
#define PARTICLESFACTORY_H

#include "Params.h"
#include "Particles.h"

class ParticlesFactory {
public:
    static Particles* create(const Params& params);
};

#endif
