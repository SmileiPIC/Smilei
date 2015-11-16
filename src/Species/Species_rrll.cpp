#include "Species_rrll.h"

#include <iostream>

#include "Particles.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrll::Species_rrll( Params& params, int ispec, Patch* patch )
    : Species( params, ispec, patch)
{
    DEBUG("Species " << ispec << "created as rrLL");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrll::~Species_rrll()
{
    DEBUG("Species rrLL deleted ");
}
