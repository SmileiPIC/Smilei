#ifndef GPU_RAND
#define GPU_RAND

#if defined(_GPU)
    #include <openacc_curand.h>
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
    #include <hiprand.hpp>
#endif

#endif