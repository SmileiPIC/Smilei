
#ifndef PROJECTOR1D_H
#define PROJECTOR1D_H

#include "Projector.h"
#include "PicParams.h"

class Projector1D : public Projector {
public:
	Projector1D(PicParams* params) : Projector(params) {};
	virtual void operator() (ElectroMagn* champs, Particle* part, double gf) = 0;
	virtual void operator() (Field* rho, Particle* part) = 0;
	
protected:
	double dx_inv_;
};

#endif

