#ifndef SPECIES_H
#define SPECIES_H

#include "Particle.h"
#include "PicParams.h"
#include "Pusher.h"
//#include "PartBoundCond.h"

#include <vector>
#include <string>
#include "PicParams.h"
#include "SmileiMPI.h"
#include "Pusher.h" 
#include "Ionization.h" 


class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class PartBoundCondMix;


//! class Species
class Species
{
 public:
	//! Species creator
	Species(PicParams*, int, SmileiMPI*);
    
	//! Species destructor
	virtual ~Species();

	//! Method returning the Particle list for the considered Species
	inline std::vector<Particle*> getParticlesList() const {return particles;}
	inline std::vector<Particle*>& getParticlesList() {return particles;}

	//! Method returning the effective number of Particles for the considered Species
	// size() = npart_effective
	inline unsigned int getNbrOfParticles() const {return particles.size();}
	// capacity() = vect ever oversize
	inline unsigned int getParticlesCapacity() const {return particles.capacity();}

	//! Method calculating the Particle dynamics (interpolation, pusher, projection)
	virtual void dynamic(double time, ElectroMagn* champs, Interpolator* interp, Projector* proj, SmileiMPI* smpi);

	//! Method used to initialize the Particle position in a given cell
	void initPosition(unsigned int, unsigned int, unsigned int *, unsigned int, std::vector<double>, std::string);

	//! Method used to initialize the Particle 3d momentum in a given cell
	void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double>&);
    
	//! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
	void initWeight(PicParams*, unsigned int, unsigned int, double);

	//! Method used to save all Particles properties for the considered Species
	void dump(std::ofstream&);

	//! Vector containing all Particles of the considered Species
	std::vector<Particle*> particles;
	//std::vector<int> index_of_particles_to_exchange;
    
	//! calculate mean charge of particles
	double meanCharge();
	
 private:
	//! species number ()
	unsigned int speciesNumber;
	
	//! Effective number of particles (different than the maximum number of particles)
	unsigned int npart_effective;

	//! Number of steps for Maxwell-Juettner cumulative function integration
	//! \todo{Put in a code constant class}
	unsigned int nE;
	
	//! Parameter used when defining the Maxwell-Juettner function (corresponds to a Maximum energy)
	double muEmax;
    
	//! Parameter used when defining the Maxwell-Juettner function (corresponds to a discretization step in energy)
	double dE;
    
	//! Number of spatial dimension for the particles
	unsigned int ndim;
	
	//! Time over which Particles of the considered Species remain frozen
	double time_frozen;

	//! Method used to apply boundary-condition for the Particles of the considered Species
	PartBoundCond* partBoundCond;
	
	//! Method used to Push the particles (change momentum & change position)
	Pusher* Push;
	
	//! Ionization method 
	Ionization* Ionize;

};

#endif
