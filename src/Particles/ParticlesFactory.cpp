// -----------------------------------------------------------------------------
//
//! \file ParticleFactory.cpp
//
//! \brief Factory to select the right Particles logic (CPU/GPU)
//
// -----------------------------------------------------------------------------
#include "ParticlesFactory.h"

#if defined( SMILEI_ACCELERATOR_GPU_OACC ) || defined( SMILEI_ACCELERATOR_GPU_OMP )
extern "C" void* CreateGPUParticles( const void* parameters,
                                     const void* a_parent_patch );
#endif

Particles* ParticlesFactory::create( const Params& parameters,
                                     const Patch&  a_parent_patch )
{
    Particles* particles = nullptr;

    if( parameters.gpu_computing ) {

        // We export a C interface to avoid potential ABI problems
        // that could occur when using two different compilers (e.g., one to
        // compile cuda/hip and another one for the host code).
#if defined( SMILEI_ACCELERATOR_GPU_OACC ) || defined( SMILEI_ACCELERATOR_GPU_OMP )
        particles = static_cast<Particles*>( CreateGPUParticles( &parameters, &a_parent_patch ) );
#else
        SMILEI_UNUSED( a_parent_patch );
        ERROR( "Unreachable state reached, Smilei was not built with GPU support!" )
#endif
    } else {
        particles = new Particles();
    }

    if( particles == nullptr ) {
        ERROR( "particles could not be allocated!" )
    }

    return particles;
}
