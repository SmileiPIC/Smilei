#include "PartCompTimeAM2Order.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

PartCompTimeAM2Order::PartCompTimeAM2Order( ) : PartCompTime() {};

float PartCompTimeAM2Order::getParticleComputationTimeVecto( const float log_particle_number ) {
    
    float r = 0;
    float x;
    
    // Skylake 8168 (Ex: Irene Joliot-Curie)
    #if defined __INTEL_SKYLAKE_8168
        r = 3.429175824327019e+00;
        x = log_particle_number;
        r +=  -1.082984887450181e+00 * x;
        x = x * log_particle_number;
        r +=  -5.343864432268729e-01 * x;
        x = x * log_particle_number;
        r += + 3.926641303143978e-01 * x;
        x = x * log_particle_number;
        r +=  -9.002576768396432e-02 * x;
        x = x * log_particle_number;
        r += + 8.923748953271823e-03 * x;
        x = x * log_particle_number;
        r +=  -3.168138264378548e-04 * x;
    // General fit
    #else
    #endif
    return r;
};

float PartCompTimeAM2Order::getParticleComputationTimeScalar( const float log_particle_number ) {
    
    float r = 0;
    float x;
    
    // Skylake 8168 (Ex: Irene)
    #if defined __INTEL_SKYLAKE_8168
        x = log_particle_number;
        r +=  -1.481007753777051e-02 * x;
    // General fit
    #else
    #endif
    
    return r;
};
