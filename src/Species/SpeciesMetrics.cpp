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
    #pragma omp simd
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
    float vecto_time_loc ( vecto_time  );
    float scalar_time_loc( scalar_time );

    double log_particle_number;
    int particle_number;

    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';

    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for (unsigned int ic=0; ic < species_loc_bmax.size(); ic++)
    {
        if (species_loc_bmax[ic] > 0)
        {
            // Max of the fit
            particle_number = std::min(species_loc_bmax[ic],256);
            // Convesion in log
            log_particle_number = log(double(particle_number));
            vecto_time_loc += SpeciesMetrics::get_particle_computation_time_vectorization(log_particle_number)*species_loc_bmax[ic];
            scalar_time_loc += SpeciesMetrics::get_particle_computation_time_scalar(log_particle_number)*species_loc_bmax[ic];
        }
    }
    vecto_time  = vecto_time_loc;
    scalar_time = scalar_time_loc;
}

// -----------------------------------------------------------------------------
//! Evaluate the time to compute all particles
//! in the current patch with vectorized operators
// -----------------------------------------------------------------------------
void SpeciesMetrics::get_computation_time(const std::vector<int> & species_loc_bmax,
                                          float & vecto_time,
                                          float & scalar_time)
{
    float vecto_time_loc ( vecto_time  );
    float scalar_time_loc( scalar_time );

    float log_particle_number;
    int particle_number;

    //std::cout << SpeciesMetrics::get_particle_computation_time_vectorization(log(32.0)) << " "
    //          << SpeciesMetrics::get_particle_computation_time_scalar(log(32.0)) << '\n';

    // Loop over the cells
    #pragma omp simd reduction(+:vecto_time_loc,scalar_time_loc) private(particle_number,log_particle_number)
    for (unsigned int ic=0; ic < species_loc_bmax.size(); ic++)
    {
        if (species_loc_bmax[ic] > 0)
        {
            // Max of the fit
            particle_number = std::min(species_loc_bmax[ic],256);
            // Convesion in log
            log_particle_number = log(float(particle_number));
            vecto_time_loc += SpeciesMetrics::get_particle_computation_time_vectorization(log_particle_number)*species_loc_bmax[ic];
            scalar_time_loc += SpeciesMetrics::get_particle_computation_time_scalar(log_particle_number)*species_loc_bmax[ic];
        }
    }

    vecto_time  = vecto_time_loc;
    scalar_time = scalar_time_loc;
}
