#ifndef SPECIES_H
#define SPECIES_H

#include "Particle.h"
#include <vector>
#include <string>
#include "PicParams.h"
#include "Pusher.h" 

//#include "PartBoundCond.h" 

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;

class Species {	
public:
	Species(PicParams*, int);
	virtual ~Species();

	inline std::vector<Particle*> getParticlesList() const {return particles;}
	inline unsigned int getNbrOfParticles() const {return npart_effective;}

	virtual void dynamic(double, ElectroMagn* , Interpolator* , Projector* );
	void initPosition(unsigned int, unsigned int, size_t *, unsigned int, std::vector<double>, std::string);
	void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double>&);
	void initMassChargeWeight(PicParams*, int, int, double);

	void dump(std::ofstream&);

private:
	std::vector<Particle*> particles;
	unsigned int npart_effective;

	//! number of steps for Maxwell-Juettner cumulative function integration
	unsigned int nE;
	
	//! maximum Energy for Maxwell-Juettner function
	double muEmax;
	double dE;
	unsigned int ndim;
	
	double time_frozen;

	PartBoundCond* partBoundCond;
	Pusher* Push;

};

#endif

