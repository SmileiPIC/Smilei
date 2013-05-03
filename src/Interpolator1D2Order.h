
#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H

#include "Interpolator1D.h"

class Interpolator1D2Order : public Interpolator1D {
public:
	Interpolator1D2Order(PicParams*, SmileiMPI*);
	void operator() (ElectroMagn* champs, Particle* part, chLocaux* ELoc, chLocaux* BLoc);
private:
};

#endif

