#include "SpeciesNorm.h"

#include <iostream>

#include "Particles.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for SpeciesNorm
// ---------------------------------------------------------------------------------------------------------------------
SpeciesNorm::SpeciesNorm( Params& params, Patch* patch )
  : Species( params, patch )
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for SpeciesNorm
// ---------------------------------------------------------------------------------------------------------------------
SpeciesNorm::~SpeciesNorm()
{
}
