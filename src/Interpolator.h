
#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "PicParams.h"
#include "SmileiMPI.h"

class ElectroMagn;
class Particle;

class Interpolator {
public:
	Interpolator(PicParams*, SmileiMPI*){;};
	virtual void operator() (ElectroMagn* champs, Particle* part, chLocaux* ELoc, chLocaux* BLoc) = 0;
private:
};

#endif

