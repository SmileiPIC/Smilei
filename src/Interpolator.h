
#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "PicParams.h"

class ElectroMagn;
class Particle;

class Interpolator {
public:
	Interpolator(PicParams *){;};
	virtual void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc) = 0;
private:
};

#endif

