
#include "Species_rrll.h"
#include "Particle.h"

#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

#include <iostream>

using namespace std;

Species_rrll::Species_rrll( PicParams* params, int ispec, SmileiMPI* smpi)
: Species( params, ispec, smpi )
{
	DEBUG(10,"Species_rrll created "<<ispec);
}

Species_rrll::~Species_rrll()
{
	DEBUG(10,"Species_rrll deleted");
}

void Species_rrll::dynamic(double time_dual, ElectroMagn* Champs, Interpolator* Interp, Projector* proj, SmileiMPI* smpi)
{
	DEBUG(10,"Species_rrll dynamic");
}

