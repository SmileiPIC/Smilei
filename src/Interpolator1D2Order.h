
#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H

#include "Interpolator1D.h"

class Interpolator1D2Order : public Interpolator1D {
public:
	Interpolator1D2Order(PicParams *);
	void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc);
private:
};

#endif

