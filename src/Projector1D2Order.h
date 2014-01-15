
#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
	Projector1D2Order(PicParams*, SmileiMPI* smpi);
	~Projector1D2Order();
	void operator() (ElectroMagn* champs, Particle* part, double gf);
	void operator() (Field* rho, Particle* part);
	void operator() (double* Jx, double* Jy, double* Jz, Particle* part, double gf, unsigned int bin, unsigned int b_dim0);
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particle* part, double gf);
private:
    double dx_ov_dt;
};

#endif

