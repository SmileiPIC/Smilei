#ifndef INTERPOLATOR1D4ORDER_H
#define INTERPOLATOR1D4ORDER_H


#include "Interpolator1D.h"
#include "Field1D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 4th order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D4Order : public Interpolator1D
{

public:
    Interpolator1D4Order(PicParams&, SmileiMPI*);
    ~Interpolator1D4Order(){};

    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc);
    inline double compute( double* coeff, Field1D* f, int idx) {
	double interp_res =  coeff[0] * (*f)(idx-2)   + coeff[1] * (*f)(idx-1)   + coeff[2] * (*f)(idx) + coeff[3] * (*f)(idx+1) + coeff[4] * (*f)(idx+2);
	return interp_res;
    };

private:
    double dble_1_ov_384 ;
    double dble_1_ov_48 ;
    double dble_1_ov_16 ;
    double dble_1_ov_12 ;
    double dble_1_ov_24 ;
    double dble_19_ov_96 ;
    double dble_11_ov_24 ;
    double dble_1_ov_4 ;
    double dble_1_ov_6 ;
    double dble_115_ov_192 ;
    double dble_5_ov_8 ;

    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Interpolation coefficient on Prim grid
    double coeffp_[5];
    // Interpolation coefficient on Dual grid
    double coeffd_[5];


};//END class

#endif
