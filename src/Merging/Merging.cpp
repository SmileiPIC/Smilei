// ----------------------------------------------------------------------------
//! \file Merging.cpp
//
//! \brief Class implementation for the generic class
//!  Merging dedicated to the particle merging.
//
//! Creation - 01/2019 - Mathieu Lobet
//
// ----------------------------------------------------------------------------

#include "Merging.h"

// -----------------------------------------------------------------------------
//! Constructor for Merging
// input: simulation parameters & Species index
//! \param params simulation parameters
//! \param species Species index
//! \param rand local random generator
// -----------------------------------------------------------------------------
Merging::Merging( Params &params, Species *species, Random * rand )
{
    // minimum particles per cell to process the merging
    min_particles_per_cell_ = species->merge_min_particles_per_cell_;
    
    // Pointer to the local patch random generator
    rand_ = rand;
}

// -----------------------------------------------------------------------------
//! Destructor for Merging
// -----------------------------------------------------------------------------
Merging::~Merging()
{
}
