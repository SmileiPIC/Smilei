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
    Interpolator2D2Order(Params&, Patch*);
    ~Interpolator2D2Order()  final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)  final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc)  final ;

    inline double compute( double* coeffx, double* coeffy, Field2D* f, int idx, int idy) {
        double interp_res(0.);
        //unroll ?
        for (int iloc=-1 ; iloc<2 ; iloc++) {
            for (int jloc=-1 ; jloc<2 ; jloc++) {
                interp_res += *(coeffx+iloc) * *(coeffy+jloc) * (*f)(idx+iloc,idy+jloc);
            }
        }
        return interp_res;
    };  

private:
    // Last prim index computed
    int ip_, jp_;
    // Last dual index computed
    int id_, jd_;
    // Last delta computed
    double deltax, deltay;
    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3];


};//END class

#endif
