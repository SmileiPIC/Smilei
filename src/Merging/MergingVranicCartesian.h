// ----------------------------------------------------------------------------
//! \file MergingVranicCartesian.h
//
//! \brief Header for the class MergingVranicCartesian
//! Particle merging with the method of Vranic et al.
//! Vranic CPC 191 65-73 (2015)
//
//! Mathieu Lobet
//! 01/2019
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#ifndef MERGINGVRANICCARTESIAN_H
#define MERGINGVRANICCARTESIAN_H

#include <cmath>

#include "Merging.h"


//------------------------------------------------------------------------------
//! MergingVranicCartesian class: holds parameters and functions to apply the
//! Vranic et al. particle merging algorithm.
//------------------------------------------------------------------------------
class MergingVranicCartesian : public Merging
{

public:

    //! Constructor for RadiationLandauLifshitz
    MergingVranicCartesian( Params &params, Species *species, Random * rand );

    //! Destructor for RadiationLandauLifshitz
    ~MergingVranicCartesian();

    // ---------------------------------------------------------------------
    //! Overloading of () operator: perform the Vranic particle merging
    //! \param particles   particle object containing the particle
    //!                    properties
    //! \param smpi        MPI properties
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param count       Final number of particles
    // ---------------------------------------------------------------------
    virtual void operator()(
        double mass,
        Particles &particles,
        std::vector <int> &mask,
        SmileiMPI *smpi,
        int istart,
        int iend,
        int & count);
        //unsigned int &remaining_particles,
        //unsigned int &merged_particles);

protected:

    // Parameters __________________________________________________

    // discretization dans chaque direction
    unsigned int dimensions_[3];

    // Minimum and maximum number of particles per packet to merge
    unsigned int min_packet_size_;
    unsigned int max_packet_size_;

    // Minimum and maximum number of particles per packet to merge
    double min_momentum_cell_length_[3];
    
    // Active the accumultion correction
    bool accumulation_correction_;

    // Discretization scale
    bool log_scale_;
    
    // Minimum momentum value in log scale
    double min_momentum_log_scale_;

private:

};

#endif
