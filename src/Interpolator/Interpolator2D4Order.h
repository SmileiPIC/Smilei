#ifndef INTERPOLATOR2D4ORDER_H
#define INTERPOLATOR2D4ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D4Order : public Interpolator2D
{

public:
    Interpolator2D4Order(Params&, Patch*);
    ~Interpolator2D4Order() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;

    inline double compute( double* coeffx, double* coeffy, Field2D* f, int idx, int idy) {
        double interp_res(0.);
        for (int iloc=-2 ; iloc<3 ; iloc++) {
            for (int jloc=-2 ; jloc<3 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * (*f)(idx+iloc,idy+jloc);
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
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[5], coeffyp_[5];
    // Interpolation coefficient on Dual grid
    double coeffxd_[5], coeffyd_[5];

};//END class

#endif
