#ifndef INTERPOLATOR2D2ORDER_H
#define INTERPOLATOR2D2ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2Order : public Interpolator2D
{

public:
    Interpolator2D2Order(PicParams&, SmileiMPI*);
    ~Interpolator2D2Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
    inline double compute( double* coeffx, double* coeffy, Field2D* f, int idx, int idy) {
	double interp_res(0.);
	for (int iloc=0 ; iloc<3 ; iloc++) {
	    for (int jloc=0 ; jloc<3 ; jloc++) {
		interp_res += coeffx[iloc] * coeffy[jloc] * (*f)(idx+iloc,idy+jloc);
	    }
	}
	return interp_res;
    };  

private:
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];


};//END class

#endif
