#ifndef PARTICLES_PARTICLESFACTORY_H
#define PARTICLES_PARTICLESFACTORY_H

#include "Params.h"
#include "Particles.h"

class ParticlesFactory {
public:
    static Particles* create(const Params& params);
};

#endif
