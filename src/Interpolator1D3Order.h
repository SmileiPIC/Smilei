#ifndef INTERPOLATOR1D3ORDER_H
#define INTERPOLATOR1D3ORDER_H

#include "Interpolator1D.h"

class Interpolator1D3Order : public Interpolator1D {
public:
    Interpolator1D3Order(PicParams *, SmileiMPI*);
    ~Interpolator1D3Order();

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
protected:
    double dble_1ov6;
    double dble_2ov3;
private:
};

#endif

