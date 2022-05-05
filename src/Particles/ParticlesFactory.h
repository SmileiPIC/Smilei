// -----------------------------------------------------------------------------
//
//! \file ParticleCreator.h
//
//! \brief Class with functions to create particles
//
// -----------------------------------------------------------------------------

#ifndef PARTICLESFACTORY_H
#define PARTICLESFACTORY_H

#include "Particles.h"
// #ifdef _GPU
// #include "nvidiaParticles.h"
// #endif
#include "Params.h"

class ParticlesFactory
{
public:
    static Particles *create( Params &params )
    {
        Particles *particles = NULL;
        
        // CPU version
        if( !params.gpu_computing ) {
            particles = new Particles();
        }
        
        // GPU version
// #ifdef _GPU
//         else {
//             particles = new nvidiaParticles();
//         }
// #endif
        return particles;
    }
};

#endif
