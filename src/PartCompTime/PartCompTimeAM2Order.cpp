#include "PartCompTimeAM2Order.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

PartCompTimeAM2Order::PartCompTimeAM2Order( ) : PartCompTime() {};

// -----------------------------------------------------------------------------
//! Evaluate the time (simple precision) to compute all particles
//! in the current patch with vectorized operators
//! @param count the numer of particles per cell
//! @aram vecto_time time in vector mode
//! @aram scalar_time time in scalar mode
// -----------------------------------------------------------------------------
void PartCompTimeAM2Order::operator()(  const std::vector<int> &count,
                                float &vecto_time,
                                float &scalar_time  )
{

    float log_particle_number;
    float particle_number;
    float vecto_time_loc = 0;
    float scalar_time_loc = 0;
        
    // std::cout << getParticleComputationTimeVecto(std::log(32.0)) << " "
    //           << getParticleComputationTimeScalar(std::log(32.0)) << '\n';
    
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
    #endif
    return r;
};

float PartCompTimeAM2Order::getParticleComputationTimeScalar( const float log_particle_number ) {
    
    float r = 0;
    
    // Skylake 8168 (Ex: Irene)
    #if defined __INTEL_SKYLAKE_8168
        r = 9.795834451761070e-01 + -1.627755115838038e-02*log_particle_number;
    // General fit
    #else
        r = 9.795834451761070e-01 + -1.627755115838038e-02*log_particle_number;
    #endif
    
    return r;
};
