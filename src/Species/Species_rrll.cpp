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
Species_rrll::Species_rrll( Params& params, SmileiMPI* smpi)
    : Species( params, smpi )
{
    particles.isRadReaction=true;
    
    ERROR("Creating a RRLL species: this is a work in progress and is still not working. Exiting");

    DEBUG("Species is being created as rrLL");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrll::~Species_rrll()
{
    DEBUG("Species rrLL deleted ");
}
