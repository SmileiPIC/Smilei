#ifndef PROJECTOR2D4ORDER_H
#define PROJECTOR2D4ORDER_H

#include "Projector2D.h"


class Projector2D4Order : public Projector2D {
public:
	Projector2D4Order(PicParams*, SmileiMPI* smpi);
	~Projector2D4Order();

	void operator() (ElectroMagn* champs, Particle* part, double gf);
	void operator() (double* Jx, double* Jy, double* Jz, Particle* part, double gf, unsigned int bin, unsigned int b_dim0);
	void operator() (Field* rho, Particle* part);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particle* part, double gf);
private:
    double one_third;
};

#endif

