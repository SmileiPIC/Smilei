
#ifndef PROJECTOR2D2ORDER_H
#define PROJECTOR2D2ORDER_H

#include "Projector2D.h"

class Projector2D2Order : public Projector2D {
public:
	Projector2D2Order(PicParams*, SmileiMPI* smpi);
	~Projector2D2Order();

	void operator() (ElectroMagn* champs, Particle* part, double gf);
	void operator() (Field* rho, Particle* part);
private:
    double dx_ov_dt;
};

#endif

