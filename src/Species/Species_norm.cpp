#include "Species_norm.h"

#include <iostream>

#include "Particles.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::Species_norm( Params& params, Patch* patch )
  : Species( params, patch )
{
    DEBUG("Species is being created as norm");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::~Species_norm()
{
    DEBUG("Species norm deleted ");
}
