
#include "Pusher.h"
#include "PicParams.h"

Pusher::Pusher(PicParams *params, int ispec)
{
	mass_   = params->species_param[ispec].mass;
	charge_ = params->species_param[ispec].charge;
	charge_over_mass_ = charge_/mass_;
	dt      = params->timestep;
	dts2    = params->timestep/2.; 

}

