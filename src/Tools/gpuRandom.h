#ifndef GPU_RANDOM
#define GPU_RANDOM

namespace smilei {

    namespace gpu {

#if defined(_GPU)
            #include <openacc_curand.h>
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
            #include <hiprand.hpp>
#endif

        class Random
        {
            public:
                
#if defined(_GPU)
            curandState_t state;
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
            hiprandState_t state;
#endif
                
        }; // end Random class definition
    
    } // end namespace gpu
} // end namespace smilei

#endif