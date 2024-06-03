#ifndef GPU_RANDOM
#define GPU_RANDOM

#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    // #include <openacc_curand.h>
    #include "curand_kernel.h"
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    // #define __HIP_PLATFORM_HCC__
    // #define __HIP_PLATFORM_AMD__
    // #include <hiprand.hpp>
    #include "Random.h"
#else
    #include "Random.h"
#endif

namespace smilei {
    namespace tools {
        namespace gpu {

            /// Universal, device/host prng.
            ///
            /// There is already a class named Random in Smilei but its not in
            /// the same namespace.
            ///
            /// Ideally we would not depend at all on openacc_curand/hiprand
            /// and use our own PRNG implementation.
            ///
            struct Random
            {
            protected:
                using State =
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                    ::curandState_t;
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                    // TODO
                    // ::hiprandState_t;
                    ::Random;
#else
                    // In src/Tools/Random.h
                    ::Random;
#endif

            public:
                Random()
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                    : a_state_{ 0xDEADBEEFU }
#else
                    : a_state_{ 0xDEADBEEFU }
#endif
                {
                    // EMPTY
                }

                // Initialization
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                void init( unsigned long long seed,
                           unsigned long long seq,
                           unsigned long long offset )
                {
                    // Cuda generator initialization
                    ::curand_init( seed, seq, offset, &a_state_ );
                }
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                void init( unsigned long long seed,
                           unsigned long long ,
                           unsigned long long  )
                {
                    // Hip generator initialization
                    // ::hiprand_init( seed, seq, offset, &state );
                    a_state_ = State{ static_cast<unsigned int>( seed ) };
                }
#else
                void init( unsigned long long seed,
                           unsigned long long ,
                           unsigned long long  )
                {
                    a_state_ = State{ static_cast<unsigned int>( seed ) };
                }
#endif

                // Initialization
                double uniform()
                {
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
                    return ::curand_uniform( &a_state_ );
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
                    // TODO
                    // return ::hiprand_uniform( &state );
                    return a_state_.uniform();
#else
                    return a_state_.uniform();
#endif
                }

                State a_state_;
            }; // end Random class definition

        } // namespace gpu
    }     // namespace tools
} // end namespace smilei

#endif
