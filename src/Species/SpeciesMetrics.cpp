/*! @file SpeciesMetrics.cpp

  @brief class SpeciesMetrics: This class contains metrics operators to evaluate
  the computation cost of a patch to treat all particles and to determine
  which type of operators should be used (vecto or not)

  @date 2018-04-20
 */

#include "SpeciesMetrics.h"



// -----------------------------------------------------------------------------
//! Return the number of cells that contain more than
//! `particle_threshold` particles
// -----------------------------------------------------------------------------
float SpeciesMetrics::get_ratio_number_of_vecto_cells( const std::vector<int> &count,
        const int particle_threshold )
{
    // - max_number_of_particles_per_cells: the maximum number of particles
    //   per cell in this patch for this species
    // - min_number_of_particles_per_cells: the minimum number of particles
    //   per cell in this patch for this species
    // Loop on all cells
    int number_of_vecto_cells = 0;
    int number_of_non_zero_cells = 0;
    float ratio_number_of_vecto_cells = 0;
    //min_number_of_particles_per_cells = count[0];
    //max_number_of_particles_per_cells = 0;
    #pragma omp simd reduction(+:number_of_vecto_cells,number_of_non_zero_cells)
    for( unsigned int ic=0; ic < count.size(); ic++ ) {
        //max_number_of_particles_per_cells = max(count[ic-1],max_number_of_particles_per_cells);
        //min_number_of_particles_per_cells = min(count[ic-1],min_number_of_particles_per_cells);
        if( count[ic-1] >= particle_threshold ) {
            number_of_vecto_cells ++;
        }
        if( count[ic-1] > 0 ) {
            number_of_non_zero_cells++;
        }
    }
    ratio_number_of_vecto_cells = float( number_of_vecto_cells ) / float( number_of_non_zero_cells );
    
    return ratio_number_of_vecto_cells;
}

// -----------------------------------------------------------------------------
//! Evaluate the time to compute all particles
//! in the current patch with vectorized operators
// -----------------------------------------------------------------------------
void SpeciesMetrics::get_computation_time( const std::vector<int> &count,
        double &vecto_time,
        double &scalar_time )
{
    double log_particle_number;
    double particle_number;
    double vecto_time_loc = 0;
    double scalar_time_loc = 0;
    
    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';
    
    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for( unsigned int ic=0; ic < count.size(); ic++ ) {
        if( count[ic] > 0 ) {
            // Max of the fit
            particle_number = std::min( double( count[ic] ), 256.0 );
            // Convesion in log
            log_particle_number = log( particle_number );
            vecto_time_loc += SpeciesMetrics::get_particle_computation_time_vectorization( log_particle_number )*count[ic];
            scalar_time_loc += SpeciesMetrics::get_particle_computation_time_scalar( log_particle_number )*count[ic];
        }
    }
    scalar_time = scalar_time_loc;
    vecto_time = vecto_time_loc;
}

// -----------------------------------------------------------------------------
//! Evaluate the time to compute all particles
//! in the current patch with vectorized operators
// -----------------------------------------------------------------------------
void SpeciesMetrics::get_computation_time( const std::vector<int> &count,
        float &vecto_time,
        float &scalar_time )
{
    float log_particle_number;
    float particle_number;
    float vecto_time_loc = 0;
    float scalar_time_loc = 0;
    
    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';
    
    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for( unsigned int ic=0; ic < count.size(); ic++ ) {
        if( count[ic] > 0 ) {
            // Max of the fit
            particle_number = std::min( float( count[ic] ), float(256.0) );
            // Convesion in log
            log_particle_number = log( particle_number );
            vecto_time_loc += get_particle_computation_time_vectorization( log_particle_number )*count[ic];
            scalar_time_loc += get_particle_computation_time_scalar( log_particle_number )*count[ic];
        }
    }
    vecto_time = vecto_time_loc;
    scalar_time = scalar_time_loc;
}

//! Evaluate the time necessary to compute `particle_number` particles
//! using vectorized operators
/*#pragma omp declare simd
double SpeciesMetrics::get_particle_computation_time_vectorization(const double log_particle_number)
{
    return 3.053068997965275e+00 + log_particle_number * (-1.346529777726451e+00 + log_particle_number * (2.262009704511124e-01 +
   log_particle_number * ( -1.220834603123080e-02 - 7.983397022180499e-05 * log_particle_number ) ) ) ;
};*/

//! Evaluate the time necessary to compute `particle_number` particles
//! using vectorized operators
//#pragma omp declare simd
float SpeciesMetrics::get_particle_computation_time_vectorization( const float log_particle_number )
{
// Cascade lake 6248 (Ex: Jean Zay)
#if defined __INTEL_CASCADELAKE_6248
    return 2.851575999756124e+00 + log_particle_number * ( -1.103410559576951e+00 + log_particle_number * ( 6.520005335065826e-02
            + log_particle_number * ( 3.143999691029673e-02 - 3.878426186471072e-03 * log_particle_number ) ) ) ;
// Skylake 8168 (Ex: Irene Joliot-Curie)
#elif defined __INTEL_SKYLAKE_8168
    return  2.873965603217334e+00 + log_particle_number * ( -1.018178658950980e+00 + log_particle_number * ( -2.390999177899332e-02
            + log_particle_number * ( 5.302690106220765e-02 - 5.500324176161280e-03 * log_particle_number ) ) ) ;
// Knight Landings Intel Xeon Phi 7250 (Ex: Frioul)
#elif defined __INTEL_KNL_7250
    return  3.391615458521049e+00 + log_particle_number * ( -1.948861281215199e+00 + log_particle_number * ( 6.609030611761257e-01
            + log_particle_number * ( -1.252595460426959e-01 + 9.287025545185804e-03 * log_particle_number ) ) ) ;
// Broadwell Intel Xeon E5-2697 v4 (Ex: Tornado)
#elif defined __INTEL_BDW_E5_2697_V4
    return  4.661824411143119e+00 + log_particle_number * ( -2.010116307618810e+00 + log_particle_number * (1.940828611778672e-01
            + log_particle_number * ( 3.249709067117774e-02 - 4.732086199743545e-03 * log_particle_number ) ) ) ;
// Haswell Intel Xeon E5-2680 v3 (Ex: Jureca)
#elif defined __INTEL_HSW_E5_2680_v3
    return  2.893485213852858e+00 + log_particle_number * ( -1.066920754145127e+00 + log_particle_number * ( 3.666171703120181e-02
            + log_particle_number * ( 3.688297004269906e-02 - 4.127980207551420e-03 * log_particle_number ) ) ) ;
// General fit
#else
    return  2.855507642982689e+00 + log_particle_number * ( -1.192593070397785e+00 + log_particle_number * ( 1.447576003168199e-01
            + log_particle_number * ( 8.410553824987992e-03 - 1.760649180606238e-03 * log_particle_number ) ) ) ;
#endif

};

//! Evaluate the time necessary to compute `particle_number` particles
//! using scalar operators
/*#pragma omp declare simd
double SpeciesMetrics::get_particle_computation_time_scalar(const double log_particle_number)
{
    return   -1.476070257489217e-02 * log_particle_number
+ 9.539747447809775e-01;
};*/

//! Evaluate the time necessary to compute `particle_number` particles
//! using scalar operators
//#pragma omp declare simd
float SpeciesMetrics::get_particle_computation_time_scalar( const float log_particle_number )
{
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
