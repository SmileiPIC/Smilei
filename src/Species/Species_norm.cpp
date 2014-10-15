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
Species_norm::Species_norm( PicParams& params, int ispec, SmileiMPI* smpi)
    : Species( params, ispec, smpi )
{
    DEBUG(20,"Species " << ispec << "created as norm");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::~Species_norm()
{
    DEBUG(20,"Species norm deleted ");
}
