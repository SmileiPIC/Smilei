
#ifndef INTERPOLATOR1D3ORDER_H
#define INTERPOLATOR1D3ORDER_H

#include "Interpolator1D.h"

class Interpolator1D3Order : public Interpolator1D {
public:
	Interpolator1D3Order(PicParams *);
	void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc);
protected:
    double dble_1ov6;
    double dble_2ov3;
private:
};

#endif

