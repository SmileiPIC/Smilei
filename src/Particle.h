#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <math.h>
#include "Tools.h"

class PicParams;

//----------------------------------------------------------------------------------------------------------------------
//! Particle class: holds the basic properties of a particle
//----------------------------------------------------------------------------------------------------------------------
class Particle {
 public:
    
	//! Constructor for Particle
	Particle(int nDim);
    
	//! Destructor for Particle
	virtual ~Particle();

	//! Method used to get the Particle position
	inline double  position( int i ) const {return Position[i];}
	//! Method used to set a new value to the Particle former position
	inline double& position( int i )       {return Position[i];}
	
	//! Method used to get the Particle position
	inline double  position_old( int i ) const {return Position_old[i];}
	//! Method used to set a new value to the Particle former position
	inline double& position_old( int i )       {return Position_old[i];}
    
	//! Method used to get the Particle momentum
	inline double  momentum( int i ) const {return Momentum[i];}
	//! Method used to set a new value to the Particle momentum
	inline double& momentum( int i )       {return Momentum[i];}
    
	//! Method used to get the Particle weight
	inline double  weight() const {return Weight[0];}
	//! Method used to set a new value to the Particle weight
	inline double& weight()       {return Weight[0];}

	//! Method used to get the Particle charge
	inline short  charge() const {return Charge[0];}
	//! Method used to set a new value to the Particle weight
	inline short& charge()       {return Charge[0];}
	
	//! Method used to get the Particle Lorentz factor
	inline double lor_fac() {return sqrt(1+pow(momentum(0),2)+pow(momentum(1),2)+pow(momentum(2),2));}
		
	//! \todo What is this doing here (MG for TV or JD)
	virtual void Print(PicParams* params);

	char *buf;
    
 private:
	//! array containing the particle position
	double* Position;
    
	//! array containing the particle former (old) positions
	double* Position_old;
    
	//! array containing the particle moments
	double*  Momentum;
    
	//! containing the particle weight: equivalent to a charge density
	double* Weight;
    
	//! charge state of the particle (multiples of e>0)
	short* Charge;
};

#endif
