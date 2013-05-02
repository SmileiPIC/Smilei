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
Species_norm::Species_norm(PicParams* params, unsigned int ispec) : Species(params, ispec)
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
