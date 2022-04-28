#ifndef GPU_RANDOM
#define GPU_RANDOM

#if defined(_GPU)
            #include <openacc_curand.h>
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
            #include <hiprand.hpp>
#endif

namespace smilei {

    namespace gpu {

        class Random
        {
            public:
            
            // State used for random generator
#if defined(_GPU)
            curandState_t state;
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
            hiprandState_t state;
#endif

            // Initialization
            template <typename T>
            static inline void init (unsigned long long seed,
                              unsigned long long seq,
                              unsigned long long offset,
                              T * state) {
#if defined(_GPU)
                curand_init(seed, seq, offset, state); //Cuda generator
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
                hiprand_init(seed, seq, offset, state); //Cuda generator initialization 
#endif
            };
                
            // Initialization
            template <typename T>
            static inline void uniform (T * state) {
#if defined(_GPU)
                curand_uniform(state); //Cuda generator
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
                hiprand_uniform(state); //Cuda generator initialization 
#endif
            };
                
        }; // end Random class definition
    
    } // end namespace gpu
} // end namespace smilei

#endif
