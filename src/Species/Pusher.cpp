#include "Pusher.h"
#include "PicParams.h"

Pusher::Pusher(PicParams& params, int ispec)
{
    mass_          = params.species_param[ispec].mass;
    one_over_mass_ = 1.0/mass_;
    dt             = params.timestep;
    dts2           = params.timestep/2.;

    nDim_          = params.nDim_particle;

}

Pusher::~Pusher()
{
}
