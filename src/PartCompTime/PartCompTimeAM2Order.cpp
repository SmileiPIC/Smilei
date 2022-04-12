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
        
    // std::cout << getParticleComputationTimeVecto(std::log(4.0)) << " "
    //           << getParticleComputationTimeScalar(std::log(4.0)) << '\n';
    
    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for( unsigned int ic=0; ic < count.size(); ic++ ) {
        if( count[ic] > 0 ) {
            // Max of the fit
            particle_number = std::min( float( count[ic] ), float(256.0) );
            // Convesion in log
            log_particle_number = std::log( particle_number );
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
    
    // Cascade lake 6248 (Ex: Jean Zay)
    #if defined __INTEL_CASCADELAKE_6248
    
        r = 3.197235443800087e+00;
        x = log_particle_number;
        r += -9.707696651705293e-01 * x;
        x = x * log_particle_number;
        r += -6.694677695526818e-01 * x;
        x = x * log_particle_number;
        r += 5.244928633089969e-01 * x;
        x = x * log_particle_number;
        r += -1.416212445676439e-01 * x;
        x = x * log_particle_number;
        r += 1.754287851325311e-02 * x;
        x = x * log_particle_number;
        r += -8.361974480680615e-04 * x;
    
    // Skylake 8168 (Ex: Irene Joliot-Curie)
    #elif defined __INTEL_SKYLAKE_8168
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
        
    // ROME 7H12 (Ex: Irene Joliot-Curie)
    #elif defined __AMD_ROME_7H12
    
        r = 4.060036149095799e+00;
        x = log_particle_number;
        r += -1.287187202618579e+00 * x;
        x = x * log_particle_number;
        r += -3.080328534429500e-01 * x;
        x = x * log_particle_number;
        r += 2.243053977668283e-01 * x;
        x = x * log_particle_number;
        r += -3.616552119492462e-02 * x;
        x = x * log_particle_number;
        r += 1.442034645317998e-03 * x;
        x = x * log_particle_number;
        r += 5.743929779052147e-05 * x;
    
    // General fit
    #else
        r = 3.562061671184216e+00;
        x = log_particle_number;
        r += -1.104209179152063e+00 * x;
        x = x * log_particle_number;
        r += -5.262022909482892e-01 * x;
        x = x * log_particle_number;
        r += 3.985374024978240e-01 * x;
        x = x * log_particle_number;
        r += -9.580833134861022e-02 * x;
        x = x * log_particle_number;
        r += 1.038743187318360e-02 * x;
        x = x * log_particle_number;
        r += -4.323215852458617e-04 * x;
    #endif
    return r;
};

float PartCompTimeAM2Order::getParticleComputationTimeScalar( const float log_particle_number ) {
    
    float r = 0;
    
    // Cascade lake 6248 (Ex: Jean Zay)
    #if defined __INTEL_CASCADELAKE_6248
    
        r = 9.733135619435868e-01 + -7.193794068651995e-03*log_particle_number;
    
    // Skylake 8168 (Ex: Irene)
    #elif defined __INTEL_SKYLAKE_8168
        r = 9.795834451761070e-01 + -1.627755115838038e-02*log_particle_number;
        
    // ROME 7H12 (Ex: Irene Joliot-Curie)
    #elif defined __AMD_ROME_7H12
        
        r = 9.082843904986886e-01 + -1.796931741585946e-02*log_particle_number;
        
    // General fit
    #else
        r = 9.537271325394608e-01 + -1.381355421429730e-02*log_particle_number;
    #endif
    
    return r;
};
