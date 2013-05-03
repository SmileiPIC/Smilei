
#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
	Projector1D2Order(PicParams*, SmileiMPI* smpi);
	void operator() (ElectroMagn* champs, Particle* part, double gf);
	void operator() (Field* rho, Particle* part);
private:
};

#endif

