
#ifndef PARTICLESFACTORY_H
#define PARTICLESFACTORY_H

#include "Particles.h"
#ifdef _GPU
#include "nvidiaParticles.h"
#endif
#include "Params.h"

class ParticlesFactory
{
public:
    static Particles *create( Params &params )
    {
        Particles *particles = NULL;
        if( !params.gpu_computing ) {
            particles = new Particles();
        }
#ifdef _GPU
        else {
            particles = new nvidiaParticles();
        }
#endif
        return particles;
    }
};

#endif
