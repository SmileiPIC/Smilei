#ifndef SPECIESMETRICS_H
#define SPECIESMETRICS_H

#include <vector>
#include <string>

//! class SpeciesMetrics: This class contains metrics operators to evaluate
//! the computation cost of a patch to treat all particles and to determine
//! which type of operators should be used (vecto or not)
class SpeciesMetrics
{

    public:

        //! Return the number of cells that contain more than `particle_threshold` particles
        static float get_ratio_number_of_vecto_cells(const std::vector<int> & species_loc_bmax,
                                                     const int particle_threshold);

    protected:



    private:


};

#endif
