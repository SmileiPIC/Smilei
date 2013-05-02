#include "Species_rrll.h"
#include "Particle.h"
#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

#include <iostream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creator for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrLL::Species_rrLL(PicParams* params, unsigned int ispec) : Species(params, ispec)
{
	DEBUG(10,"Species " << ispec << "created as rrLL");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrLL::~Species_rrLL()
{
	DEBUG(10,"Species rrLL deleted ");
}
