
#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H

#include "Interpolator.h"

class Interpolator1D : public Interpolator {
public:
	Interpolator1D(PicParams * params): Interpolator(params){;};
	virtual void operator() (ElectroMagn* champs, Particle* part, chLocaux* ELoc, chLocaux* BLoc) = 0;
protected:
	double dx_inv_;
};

#endif

