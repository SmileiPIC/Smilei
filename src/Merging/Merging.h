// ----------------------------------------------------------------------------
//! \file Merging.h
//
//! \brief Header for the generic class Merging
//! dedicated to the particle merging.
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#ifndef MERGING_H
#define MERGING_H

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "Random.h"

//  ----------------------------------------------------------------------------
//! Class Merging
//  ----------------------------------------------------------------------------
class Merging
{
public:

    //! Creator for Merging
    Merging( Params &params, Species *species, Random * rand );

    virtual ~Merging();

    //! Overloading of () operator
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    virtual void operator()(
        double mass,
        Particles &particles,
        std::vector <int> &mask,
        SmileiMPI *smpi,
        int istart,
        int iend,
        int & count) = 0;

    // parameters _______________________________________________

protected:
    
    // Local rand generator
    Random * rand_;
    
    // Minimum number of particles per cell to process the merging
    unsigned int min_particles_per_cell_;
    
private:
    
};

#endif
