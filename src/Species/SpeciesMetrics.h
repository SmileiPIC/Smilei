#ifndef SPECIESMETRICS_H
#define SPECIESMETRICS_H

#include <cstdlib>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

//! class SpeciesMetrics: This class contains metrics operators to evaluate
//! the computation cost of a patch to treat all particles and to determine
//! which type of operators should be used (vecto or not)
class SpeciesMetrics
{

    public:

        //! Return the number of cells that contain more than `particle_threshold` particles
        static float get_ratio_number_of_vecto_cells(const std::vector<int> & species_loc_bmax,
                                                     const int particle_threshold);

        //! Evaluate the time to compute all particles in the current patch with vectorized operators
        static void get_computation_time(const std::vector<int> & species_loc_bmax,
                                         double & vecto_time,
                                         double & scalar_time);

    protected:

        //! Evaluate the time necessary to compute `particle_number` particles
        //! using vectorized operators
        static double inline get_particle_computation_time_vectorization(const double log_particle_number)
        {
            return -7.898122120999234e-03 * pow(log_particle_number,4)
                    + 7.194816050066996e-02 * pow(log_particle_number,3)
                    -3.932507419636225e-02 * pow(log_particle_number,2)
                    -1.202283684111388 * log_particle_number
                    + 3.338344755724750;
        };

        //! Evaluate the time necessary to compute `particle_number` particles
        //! using scalar operators
        static double inline get_particle_computation_time_scalar(const double log_particle_number)
        {
            return -2.614094907002057e-02 * log_particle_number
                    + 9.776048513444786e-01;
        };

    private:


};

#endif
