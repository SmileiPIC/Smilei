
#include "Species_norm.h"
#include "Particle.h"

#include "Interpolator.h"
#include "Projector.h"
#include "Pusher.h"

#include <iostream>

using namespace std;

Species_norm::Species_norm( PicParams* params, int ispec, SmileiMPI* smpi)
 : Species( params, ispec, smpi )
{
	DEBUG(10,"Species norm created "<<ispec);
}

Species_norm::~Species_norm()
{
	DEBUG(10,"Species norm deleted ");
}

/*void Species_norm::dynamic(ElectroMagn* Champs, Pusher* ppush, Interpolator* Interp, Projector* proj)
{
	DEBUG(10,"dynamic Species_norm");

}*/

