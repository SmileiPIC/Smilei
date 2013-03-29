
#ifndef ParticleRad_H
#define ParticleRad_H

#include "Particle.h"

class PicParams;

class ParticleRad : public Particle {
public:	
	ParticleRad(int);
	~ParticleRad();

	inline double  radPower() const {return rad_power;};
	inline double& radPower()       {return rad_power;};
	inline double  omegaCrit() const {return omega_crit;};
	inline double& omegaCrit()       {return omega_crit;};
	
	//void Initialize(PicParams* params) {;};  
	//void Print(PicParams* params) {;};
  
private:
	double rad_power;
	double omega_crit;

};

//ostream& operator <<(std::ostream& Stream, const Particle1D& Obj)
//{
//    Stream << "sizeof Particle : " << sizeof( Obj ) << endl;
//    //Stream <<  &Obj.Psm_[0] << " " << &Obj.weight_ << " " << &Obj.masse_ << " " << &Obj.charge_;
//    Stream <<  &Obj.Psm_[0] - &Obj.weight_ << " " << &Obj.weight_ - &Obj.masse_ << " " << &Obj.masse_ - &Obj.charge_ << &Obj.charge_ - &Obj.pos_[0];
//    return Stream; // N'oubliez pas de renvoyer le flux, afin de pouvoir chaîner les appels
//}

#endif

