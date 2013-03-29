
#ifndef PROJECTOR2D_H
#define PROJECTOR2D_H

#include "Projector.h"

class Projector2D : public Projector {
public:
	virtual void operator() (ElectroMagn* champs, Particle* part, double gf) = 0;
	virtual void operator() (Field* rho, Particle* part) = 0;

private:
};

#endif

