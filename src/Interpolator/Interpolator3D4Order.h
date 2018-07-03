#ifndef INTERPOLATOR3D4ORDER_H
#define INTERPOLATOR3D4ORDER_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D4Order : public Interpolator3D
{

public:
    Interpolator3D4Order(Params&, Patch*);
    ~Interpolator3D4Order() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;

    inline double compute( double* coeffx, double* coeffy, double* coeffz, Field3D* f, int idx, int idy, int idz) {
	double interp_res(0.);
        //unroll ?
	for (int iloc=-2 ; iloc<3 ; iloc++) {
	    for (int jloc=-2 ; jloc<3 ; jloc++) {
                for (int kloc=-2 ; kloc<3 ; kloc++) {
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * *(coeffz+kloc) * (*f)(idx+iloc,idy+jloc,idz+kloc);
                }
	    }
	}
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
    int ip_, jp_, kp_;
    // Last dual index computed
    int id_, jd_, kd_;
    // Last delta computed
    double deltax, deltay, deltaz;
    // Interpolation coefficient on Prim grid
    double coeffxp_[5], coeffyp_[5], coeffzp_[5];
    // Interpolation coefficient on Dual grid
    double coeffxd_[5], coeffyd_[5], coeffzd_[5];


};//END class

#endif
