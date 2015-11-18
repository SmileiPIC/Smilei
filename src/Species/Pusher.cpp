#include "Pusher.h"
#include "Params.h"
#include "Species.h"

Pusher::Pusher(Params& params, Species *species)
{
    mass_          = species->mass;
    one_over_mass_ = 1.0/mass_;
    dt             = params.timestep;
    dts2           = params.timestep/2.;

    nDim_          = params.nDim_particle;

}

Pusher::~Pusher()
{
}
