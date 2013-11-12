
#ifndef PARTICLEFACTORY_H
#define PARTICLEFACTORY_H

#include "Particle.h"
#include "ParticleRad.h"

#include "PicParams.h"

#include "Tools.h"

class ParticleFactory {
public:
	static Particle* create(PicParams* params, int ispec) {
		Particle* part = NULL;
	    if (params->species_param[ispec].radiating) {
		    part=new ParticleRad(params->nDim_particle);
	    }
	    else {
		    part=new Particle(params->nDim_particle);
	    }

		return part;
	}

	static std::vector<Particle*> createVector(PicParams* params, int ispec, int npart_effective) {
		std::vector<Particle*> vecParticles;
		vecParticles.resize(npart_effective);

		for (int ipart=0 ; ipart<npart_effective ; ipart++) {
			vecParticles[ipart] = ParticleFactory::create(params, ispec);

		} // END for ipart

		return vecParticles;
	}

};

#endif
