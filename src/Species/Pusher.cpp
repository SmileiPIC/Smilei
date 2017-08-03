#include "Pusher.h"
#include "Params.h"
#include "Species.h"

Pusher::Pusher(Params& params, Species *species)
{
    mass_          = species->mass;
    if (mass_ > 0.)
    {
        one_over_mass_ = 1.0/mass_;
    }
    else
    {
        one_over_mass_ = 0.;
    }
    dt             = params.timestep;
    dts2           = params.timestep/2.;

    nDim_          = params.nDim_particle;

}

Pusher::~Pusher()
{
}
