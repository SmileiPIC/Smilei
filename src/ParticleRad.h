#ifndef ParticleRad_H
#define ParticleRad_H

#include "Particle.h"


class PicParams;



//----------------------------------------------------------------------------------------------------------------------
//! Class ParticleRad: contains properties of a radiating Particle
//----------------------------------------------------------------------------------------------------------------------
class ParticleRad : public Particle {
public:
    
    //! Constructor for ParticleRad
	ParticleRad(int);
    
    //! Destructor for ParticleRad
	~ParticleRad();

    //! Method used to get the power radiated by the particle
	inline double  radPower() const {return rad_power;};

    //! Method used to set a new value to the power radiated by the particle
	inline double& radPower()       {return rad_power;};

	//! Method used to get the critical frequency of the radiation emitted by the particle
    inline double  omegaCrit() const {return omega_crit;};
    
    //! Method used to set a new value to the critical frequency of the radiation emitted by the particle
	inline double& omegaCrit()       {return omega_crit;};
    
  
private:
    //! \todo For radiating particle use the eta parameter (MG)
    
    //! Power radiated away by the particle
	double rad_power;
    
    //! Critical frequency of the radiation emitted by the particle
	double omega_crit;

};


#endif
