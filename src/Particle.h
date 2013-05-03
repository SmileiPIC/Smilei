
#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include "Tools.h"

class PicParams;

//! this holds the basic properties of a particle
class Particle {
public:	
	Particle(int nDim);
	virtual ~Particle();

	inline double  position( int i ) const {return Pos_[i];}
	inline double& position( int i )       {return Pos_[i];}
	inline double  moments( int i ) const {return Psm_[i];}
	inline double& moments( int i )       {return Psm_[i];}
	inline double  chargeDensity() const {return charge_density_[0];}
	inline double& chargeDensity()       {return charge_density_[0];}
	
	//virtual void Initialize(PicParams*, int, int, int, int, int, double, double, double, double, double, double, double);
	virtual void Print(PicParams* params);

	//void* operator new(size_t);
	//void operator delete(void*);
	double *buf;//  = new char[sizeof(Particle)];

	
private:
	//! array of positions
	double* Pos_;
	//! array of moments
	double*  Psm_;
	//! weight*charge \todo{specify units!}
	double* charge_density_;

};

#endif

