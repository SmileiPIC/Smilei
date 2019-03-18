// ----------------------------------------------------------------------------
//! \file MergingVranic.h
//
//! \brief Header for the class MergingVranic
//! Particle merging with the method of Vranic et al.
//! Vranic CPC 191 65-73 (2015)
//
//! Mathieu Lobet
//! 01/2019
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#ifndef MERGINGVRANIC_H
#define MERGINGVRANIC_H

#include <cmath>

#include "Merging.h"


//------------------------------------------------------------------------------
//! MergingVranic class: holds parameters and functions to apply the
//! Vranic et al. particle merging algorithm.
//------------------------------------------------------------------------------
class MergingVranic : public Merging
{

public:

    //! Constructor for RadiationLandauLifshitz
    MergingVranic( Params &params, Species *species );

    //! Destructor for RadiationLandauLifshitz
    ~MergingVranic();

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

    // Parameters __________________________________________________

    // discretization dans chaque direction
    unsigned int dimensions_[3];

    // Minimum and maximum number of particles per packet to merge
    unsigned int min_packet_size_;
    unsigned int max_packet_size_;

protected:


private:

};

#endif
