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

#include "Merging.h"

//------------------------------------------------------------------------------
//! MergingVranic class: holds parameters and functions to apply the
//! Vranic et al. particle merging algorithm.
//------------------------------------------------------------------------------
class MergingVranic : public Merging {

public:

    //! Constructor for RadiationLandauLifshitz
    MergingVranic(Params& params, Species * species);

    //! Destructor for RadiationLandauLifshitz
    ~MergingVranic();

    // ---------------------------------------------------------------------
    //! Overloading of () operator: perform the Vranic particle merging
    //! \param particles   particle object containing the particle
    //!                    properties
    //! \param smpi        MPI properties
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    // ---------------------------------------------------------------------
    virtual void operator() (
            Particles &particles,
            SmileiMPI* smpi,
            int istart,
            int iend,
            int ithread,
            int ipart_ref = 0);

protected:


private:

};

#endif
