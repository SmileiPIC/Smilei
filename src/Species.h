#ifndef SPECIES_H
#define SPECIES_H

#include "Particle.h"
#include <vector>
#include <string>
#include "PicParams.h"
#include "SmileiMPI.h"
#include "Pusher.h" 

//#include "PartBoundCond.h" 

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;

class Species {	
public:
	Species(PicParams*, int, SmileiMPI*);
	virtual ~Species();

	inline std::vector<Particle*> getParticlesList() const {return particles;}
	inline std::vector<Particle*>& getParticlesList() {return particles;}
	// size() = npart_effective
	inline unsigned int getNbrOfParticles() const {return particles.size();}
	// capacity() = vect ever oversize
	inline unsigned int getParticlesCapacity() const {return particles.capacity();}

	void initPosition(unsigned int, unsigned int, size_t *, unsigned int, std::vector<double>, std::string);
	void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double>&);
	void initMassChargeWeight(PicParams*, int, int, double);

	virtual void dynamic(double time, ElectroMagn* champs, Interpolator* interp, Projector* proj, SmileiMPI* smpi);


	void dump(std::ofstream&);
	std::vector<Particle*> particles;

private:
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

