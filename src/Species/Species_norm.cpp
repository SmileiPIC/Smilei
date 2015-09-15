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
Species_norm::Species_norm( Params& params,  SpeciesStructure& sparams, SmileiMPI* smpi)
    : Species( params, sparams, smpi )
{
    DEBUG("Species " << sparams.species_type << "created as norm");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::~Species_norm()
{
    DEBUG("Species norm deleted ");
}
