#ifndef GPU_RAND
#define GPU_RAND

#if defined(_GPU)
    #include <hiprand.hpp>
#elif defined(SMILEI_ACCELERATOR_GPU_OMP)
    #include <openacc_curand.h>
#endif

#endif