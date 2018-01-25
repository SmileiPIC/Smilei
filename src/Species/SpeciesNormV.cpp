#include "SpeciesNormV.h"

#include <iostream>

#include "Particles.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for SpeciesNorm
// ---------------------------------------------------------------------------------------------------------------------
SpeciesNormV::SpeciesNormV( Params& params, Patch* patch )
  : SpeciesV( params, patch )
{
    DEBUG("Species is being created as norm");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for SpeciesNorm
// ---------------------------------------------------------------------------------------------------------------------
SpeciesNormV::~SpeciesNormV()
{
    DEBUG("Species norm deleted ");
}
