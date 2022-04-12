#include "PartCompTime2D2Order.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

PartCompTime2D2Order::PartCompTime2D2Order( ) : PartCompTime() {};

// -----------------------------------------------------------------------------
//! Evaluate the time (simple precision) to compute all particles
//! in the current patch with vectorized operators
//! @param count the numer of particles per cell
//! @aram vecto_time time in vector mode
//! @aram scalar_time time in scalar mode
// -----------------------------------------------------------------------------
void PartCompTime2D2Order::operator()(  const std::vector<int> &count,
                                float &vecto_time,
                                float &scalar_time  )
{

    float log_particle_number;
    float particle_number;
    float vecto_time_loc = 0;
    float scalar_time_loc = 0;
    
    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for( unsigned int ic=0; ic < count.size(); ic++ ) {
        if( count[ic] > 0 ) {
            // Max of the fit
            particle_number = std::min( float( count[ic] ), float(256.0) );
            // Convesion in log
            log_particle_number = log( particle_number );
            vecto_time_loc += getParticleComputationTimeVecto( log_particle_number )*count[ic];
            scalar_time_loc += getParticleComputationTimeScalar( log_particle_number )*count[ic];
        }
    }
    scalar_time = scalar_time_loc;
    vecto_time = vecto_time_loc;
};

float PartCompTime2D2Order::getParticleComputationTimeVecto( const float log_particle_number ) {
    
    float r = 0;
    float x;
    
    // Skylake 8168 (Ex: Irene Joliot-Curie)
    #if defined __INTEL_SKYLAKE_8168
        r = 3.763654995524249e+00;
        x = log_particle_number;
        r += -2.757580191504306e+00 * x;
        x = x * log_particle_number;
        r += 8.514971254090451e-01 * x;
        x = x * log_particle_number;
        r += -1.037213025203770e-01 * x;
        x = x * log_particle_number;
        r += 2.988386432602295e-03 * x;
        x = x * log_particle_number;
        r += 1.746776707365379e-04 * x;
        
    // ROME 7H12 (Ex: Irene Joliot-Curie)
    #elif defined __AMD_ROME_7H12
        
    // General fit
    #else
        r = 3.763654995524249e+00;
        x = log_particle_number;
        r += -2.757580191504306e+00 * x;
        x = x * log_particle_number;
        r += 8.514971254090451e-01 * x;
        x = x * log_particle_number;
        r += -1.037213025203770e-01 * x;
        x = x * log_particle_number;
        r += 2.988386432602295e-03 * x;
        x = x * log_particle_number;
        r += 1.746776707365379e-04 * x;
    #endif
    return r;
};

float PartCompTime2D2Order::getParticleComputationTimeScalar( const float log_particle_number ) {
    
    float r = 0;
    
    // Skylake 8168 (Ex: Irene)
    #if defined __INTEL_SKYLAKE_8168
        r = 9.449407213258362e-01 + 1.416491664909733e-02*log_particle_number;
    // ROME 7H12 (Ex: Irene Joliot-Curie)
    #elif defined __AMD_ROME_7H12
    // General fit
    #else
        r = 9.449407213258362e-01 + 1.416491664909733e-02*log_particle_number;
    #endif
    
    return r;
};
