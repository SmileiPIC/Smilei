#include "ParticleRad.h"

#include <iostream>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ParticleRad (used for radiating Particles)
// ---------------------------------------------------------------------------------------------------------------------
ParticleRad::ParticleRad(int nDim) : Particle(nDim)
{
	rad_power=0.0;
	omega_crit=0.0;
	DEBUG(10, "Particle rad created " << nDim << "D");
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ParticleRad
// ---------------------------------------------------------------------------------------------------------------------
ParticleRad::~ParticleRad()
{
}
