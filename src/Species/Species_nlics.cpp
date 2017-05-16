// ---------------------------------------------------------------------------------------------------------------------
//
//! \file Species_nlics.cpp
//
//! \brief Species_nlics.cpp  generic class for the species that use the Monte-Carlo 
//!  pusher for the Non Linear Inverse Compton Scattering.
//
//! \date 2017-05-12
//
// ---------------------------------------------------------------------------------------------------------------------

#include "Species_nlics.h"

#include <iostream>

#include "Particles.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species_nlics
// ---------------------------------------------------------------------------------------------------------------------
Species_nlics::Species_nlics( Params& params, Patch* patch, NLICompton &)
    : Species( params, patch )
{
    particles->isRadReaction=true;
    
    DEBUG("Species is being created as nlics");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_nlics::~Species_nlics()
{
    DEBUG("Species nlics deleted ");
}
