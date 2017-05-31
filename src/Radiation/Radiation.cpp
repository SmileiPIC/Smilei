// ----------------------------------------------------------------------------
//! \file Radiation.cpp
//
//! \brief This file contains the class functions for the generic class
//  Radiation for the particle radiation losses.
//
// ----------------------------------------------------------------------------

#include "Radiation.h"

// ---------------------------------------------------------------------------------------------------------------------
//! Constructor for Radiation
// input: simulation parameters & Species index
//! \param params simulation parameters
//! \param species Species index
// ---------------------------------------------------------------------------------------------------------------------
Radiation::Radiation(Params& params, Species * species)
{
    // Dimension position
    nDim_ = params.nDim_particle;

    // Time step
    dt    = params.timestep;

    // Inverse of the species mass
    one_over_mass_ = 1./species->mass;

    // Normalized Schwinger Electric Field
    norm_E_Schwinger = electron_mass*c_vacuum*c_vacuum
                     / (red_planck_cst*params.referenceAngularFrequency_SI);

    // Inverse of norm_E_Schwinger
    inv_norm_E_Schwinger = 1./norm_E_Schwinger;
}

// ---------------------------------------------------------------------------------------------------------------------
//! Destrructor for Radiation
// ---------------------------------------------------------------------------------------------------------------------
Radiation::~Radiation()
{
}
