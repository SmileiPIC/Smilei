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
Species_rrll::Species_rrll( PicParams& params, int ispec, SmileiMPI* smpi)
    : Species( params, ispec, smpi )
{
    DEBUG(20,"Species " << ispec << "created as rrLL");
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Species_rrLL
// ---------------------------------------------------------------------------------------------------------------------
Species_rrll::~Species_rrll()
{
    DEBUG(20,"Species rrLL deleted ");
}
