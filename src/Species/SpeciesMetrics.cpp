/*! @file SpeciesMetrics.cpp

  @brief class SpeciesMetrics: This class contains metrics operators to evaluate
  the computation cost of a patch to treat all particles and to determine
  which type of operators should be used (vecto or not)

  @date 2018-04-20
 */

#include "SpeciesMetrics.h"

//! Return the number of cells that contain more than `particle_threshold` particles
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
    for (unsigned int ic=1; ic < species_loc_bmax.size(); ic++)
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
