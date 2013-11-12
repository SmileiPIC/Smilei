#include "Species_norm.h"
#include "Particle.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

#include <iostream>

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::Species_norm( PicParams* params, int ispec, SmileiMPI* smpi)
 : Species( params, ispec, smpi )
{
	DEBUG(10,"Species " << ispec << "created as norm");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_norm
// ---------------------------------------------------------------------------------------------------------------------
Species_norm::~Species_norm()
{
	DEBUG(10,"Species norm deleted ");
}
