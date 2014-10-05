#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"

#include "Field1D.h"
//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order : public Interpolator1D
{

public:
    Interpolator1D2Order(PicParams&, SmileiMPI*);
    ~Interpolator1D2Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);

    inline double compute( double* coeff, Field1D* f, int idx) {
	double interp_res =  coeff[0] * (*f)(idx-1)   + coeff[1] * (*f)(idx)   + coeff[2] * (*f)(idx+1);
	return interp_res;
    };

private:
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Interpolation coefficient on Prim grid
    double coeffp_[3];
    // Interpolation coefficient on Dual grid
    double coeffd_[3];


};//END class

#endif
