#ifndef INTERPOLATOR1D3ORDER_H
#define INTERPOLATOR1D3ORDER_H

#include "Interpolator1D.h"
#include "Field1D.h"

class Interpolator1D3Order : public Interpolator1D {
public:
    Interpolator1D3Order(PicParams&, SmileiMPI*);
    ~Interpolator1D3Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
    inline double compute( double* coeff, Field1D* f, int idx) {
	double interp_res =  coeff[0] * (*f)(idx-1)   + coeff[1] * (*f)(idx)   + coeff[2] * (*f)(idx+1) + coeff[3] * (*f)(idx+2);
	return interp_res;
    };

protected:
    double dble_1ov6;
    double dble_2ov3;
private:
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Interpolation coefficient on Prim grid
    double coeffp_[4];
    // Interpolation coefficient on Dual grid
    double coeffd_[4];

};

#endif

