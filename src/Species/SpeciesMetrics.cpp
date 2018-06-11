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
float SpeciesMetrics::get_ratio_number_of_vecto_cells(const std::vector<int> & species_loc_bmax,
                                              const int particle_threshold)
{
    // - max_number_of_particles_per_cells: the maximum number of particles
    //   per cell in this patch for this species
    // - min_number_of_particles_per_cells: the minimum number of particles
    //   per cell in this patch for this species
    // Loop on all cells
    int number_of_vecto_cells = 0;
    int number_of_non_zero_cells = 0;
    float ratio_number_of_vecto_cells = 0;
    //min_number_of_particles_per_cells = species_loc_bmax[0];
    //max_number_of_particles_per_cells = 0;
    #pragma omp simd reduction(+:number_of_vecto_cells,number_of_non_zero_cells)
    for (unsigned int ic=0; ic < species_loc_bmax.size(); ic++)
    {
        //max_number_of_particles_per_cells = max(species_loc_bmax[ic-1],max_number_of_particles_per_cells);
        //min_number_of_particles_per_cells = min(species_loc_bmax[ic-1],min_number_of_particles_per_cells);
        if (species_loc_bmax[ic-1] >= particle_threshold)
        {
            number_of_vecto_cells ++;
        }
        if (species_loc_bmax[ic-1] > 0)
        {
            number_of_non_zero_cells++;
        }
    }
    ratio_number_of_vecto_cells = float(number_of_vecto_cells) / float(number_of_non_zero_cells);

    return ratio_number_of_vecto_cells;
}

// -----------------------------------------------------------------------------
//! Evaluate the time to compute all particles
//! in the current patch with vectorized operators
// -----------------------------------------------------------------------------
void SpeciesMetrics::get_computation_time(const std::vector<int> & species_loc_bmax,
                                          double & vecto_time,
                                          double & scalar_time)
{
    double log_particle_number;
    double particle_number;
    double vecto_time_loc = 0;
    double scalar_time_loc = 0;

    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';

    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for (unsigned int ic=0; ic < species_loc_bmax.size(); ic++)
    {
        if (species_loc_bmax[ic] > 0)
        {
            // Max of the fit
            particle_number = fmin(double(species_loc_bmax[ic]),256.0);
            // Convesion in log
            log_particle_number = log(particle_number);
            vecto_time_loc += SpeciesMetrics::get_particle_computation_time_vectorization(log_particle_number)*species_loc_bmax[ic];
            scalar_time_loc += SpeciesMetrics::get_particle_computation_time_scalar(log_particle_number)*species_loc_bmax[ic];
        }
    }
    scalar_time = scalar_time_loc;
    vecto_time = vecto_time_loc;
}

// -----------------------------------------------------------------------------
//! Evaluate the time to compute all particles
//! in the current patch with vectorized operators
// -----------------------------------------------------------------------------
void SpeciesMetrics::get_computation_time(const std::vector<int> & species_loc_bmax,
                                          float & vecto_time,
                                          float & scalar_time)
{
    float log_particle_number;
    float particle_number;
    float vecto_time_loc = 0;
    float scalar_time_loc = 0;

    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';

    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for (unsigned int ic=0; ic < species_loc_bmax.size(); ic++)
    {
        if (species_loc_bmax[ic] > 0)
        {
            // Max of the fit
            particle_number = fmin(float(species_loc_bmax[ic]),256.0);
            // Convesion in log
            log_particle_number = log(particle_number);
            vecto_time_loc += get_particle_computation_time_vectorization(log_particle_number)*species_loc_bmax[ic];
            scalar_time_loc += get_particle_computation_time_scalar(log_particle_number)*species_loc_bmax[ic];
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
    return   -7.983397022180499e-05 * pow(log_particle_number,4)
 -1.220834603123080e-02 * pow(log_particle_number,3)
+ 2.262009704511124e-01 * pow(log_particle_number,2)
 -1.346529777726451e+00 * log_particle_number
+ 3.053068997965275e+00;
};*/

//! Evaluate the time necessary to compute `particle_number` particles
//! using vectorized operators
#pragma omp declare simd
float SpeciesMetrics::get_particle_computation_time_vectorization(const float log_particle_number)
{
// Skylake 8168 (Ex: Irene)
#if defined INTEL_SKYLAKE_8168
    return    -5.500324176161280e-03 * pow(log_particle_number,4)
+ 5.302690106220765e-02 * pow(log_particle_number,3)
 -2.390999177899332e-02 * pow(log_particle_number,2)
 -1.018178658950980e+00 * log_particle_number
+ 2.873965603217334e+00;
// Knight Landings Intel Xeon Phi 7250 (Ex: Frioul)
#elif defined INTEL_KNL_7250
    return    + 9.287025545185804e-03 * pow(log_particle_number,4)
 -1.252595460426959e-01 * pow(log_particle_number,3)
+ 6.609030611761257e-01 * pow(log_particle_number,2)
 -1.948861281215199e+00 * log_particle_number
+ 3.391615458521049e+00
// General fit
#else
    return   -7.983397022180499e-05 * pow(log_particle_number,4)
 -1.220834603123080e-02 * pow(log_particle_number,3)
+ 2.262009704511124e-01 * pow(log_particle_number,2)
 -1.346529777726451e+00 * log_particle_number
+ 3.053068997965275e+00;
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
#pragma omp declare simd
float SpeciesMetrics::get_particle_computation_time_scalar(const float log_particle_number)
{
// Skylake 8168 (Ex: Irene)
#if defined INTEL_SKYLAKE_8168
    return   -1.476070257489217e-02 * log_particle_number
+ 9.539747447809775e-01;
// Knight Landings Intel Xeon Phi 7250 (Ex: Frioul)
#elif defined INTEL_KNL_7250
    return   -1.693420314189753e-02 * log_particle_number
+ 9.640406193625433e-01
// General fit
#else
    return   -3.227685432492503e-02 * log_particle_number
+ 9.344604887689714e-01;
#endif
};
