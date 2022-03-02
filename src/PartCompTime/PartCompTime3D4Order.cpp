#include "PartCompTime3D4Order.h"

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

PartCompTime3D4Order::PartCompTime3D4Order( ) : PartCompTime() {};

// -----------------------------------------------------------------------------
//! Evaluate the time (simple precision) to compute all particles
//! in the current patch with vectorized operators
//! @param count the numer of particles per cell
//! @aram vecto_time time in vector mode
//! @aram scalar_time time in scalar mode
// -----------------------------------------------------------------------------
void PartCompTime3D4Order::operator()(  const std::vector<int> &count,
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

float PartCompTime3D4Order::getParticleComputationTimeVecto( const float log_particle_number ) {
    // Cascade lake 6248 (Ex: Jean Zay)
    #if defined __INTEL_CASCADELAKE_6248
        return -3.878426186471072e-03 * pow(log_particle_number,4)
                + 3.143999691029673e-02 * pow(log_particle_number,3)
                + 6.520005335065826e-02 * pow(log_particle_number,2)
                -1.103410559576951e+00 * log_particle_number
                + 2.851575999756124e+00;
    // Skylake 8168 (Ex: Irene Joliot-Curie)
    #elif defined __INTEL_SKYLAKE_8168
        return    -5.500324176161280e-03 * pow( log_particle_number, 4 )
                  + 5.302690106220765e-02 * pow( log_particle_number, 3 )
                  -2.390999177899332e-02 * pow( log_particle_number, 2 )
                  -1.018178658950980e+00 * log_particle_number
                  + 2.873965603217334e+00;
    // Knight Landings Intel Xeon Phi 7250 (Ex: Frioul)
    #elif defined __INTEL_KNL_7250
        return    + 9.287025545185804e-03 * pow( log_particle_number, 4 )
                  -1.252595460426959e-01 * pow( log_particle_number, 3 )
                  + 6.609030611761257e-01 * pow( log_particle_number, 2 )
                  -1.948861281215199e+00 * log_particle_number
                  + 3.391615458521049e+00;
    // Broadwell Intel Xeon E5-2697 v4 (Ex: Tornado)
    #elif defined __INTEL_BDW_E5_2697_V4
        return     -4.732086199743545e-03 * pow( log_particle_number, 4 )
                   + 3.249709067117774e-02 * pow( log_particle_number, 3 )
                   + 1.940828611778672e-01 * pow( log_particle_number, 2 )
                   -2.010116307618810e+00 * log_particle_number
                   + 4.661824411143119e+00;
    // Haswell Intel Xeon E5-2680 v3 (Ex: Jureca)
    #elif defined __INTEL_HSW_E5_2680_v3
        return     -4.127980207551420e-03 * pow( log_particle_number, 4 )
                   + 3.688297004269906e-02 * pow( log_particle_number, 3 )
                   + 3.666171703120181e-02 * pow( log_particle_number, 2 )
                   -1.066920754145127e+00 * log_particle_number
                   + 2.893485213852858e+00;
    // General fit
    #else
        return    -1.760649180606238e-03 * pow(log_particle_number,4)
                + 8.410553824987992e-03 * pow(log_particle_number,3)
                + 1.447576003168199e-01 * pow(log_particle_number,2)
                 -1.192593070397785e+00 * log_particle_number
                + 2.855507642982689e+00;
    #endif
};

float PartCompTime3D4Order::getParticleComputationTimeScalar( const float log_particle_number ) {
    // Cascade lake 6248 (Ex: Jean Zay)
    #if defined __INTEL_CASCADELAKE_6248
        return  -1.109419407609368e-02 * log_particle_number
                + 9.564834436679180e-01;
    // Skylake 8168 (Ex: Irene)
    #elif defined __INTEL_SKYLAKE_8168
        return   -1.476070257489217e-02 * log_particle_number
                 + 9.539747447809775e-01;
    // Knight Landings Intel Xeon Phi 7250 (Ex: Frioul)
    #elif defined __INTEL_KNL_7250
        return   -1.693420314189753e-02 * log_particle_number
                 + 9.640406193625433e-01;
    // Broadwell Intel Xeon E5-2697 v4 (Ex: Tornado)
    #elif defined __INTEL_BDW_E5_2697_V4
        return   + 6.694852027937652e-03 * log_particle_number
                 + 9.382353109818060e-01;
    // Haswell Intel Xeon E5-2680 v3 (Ex: Jureca)
    #elif defined __INTEL_HSW_E5_2680_v3
        return    -1.716273243387051e-02 * log_particle_number
                  + 9.761935025470106e-01;
    // General fit
    #else
        return    -8.274579739804833e-03 * log_particle_number
                + 9.405673399529412e-01;
    #endif
};
