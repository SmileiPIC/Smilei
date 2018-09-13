#ifndef INTERPOLATOR2D2ORDER_H
#define INTERPOLATOR2D2ORDER_H


#include "Interpolator2D.h"
#include "Field2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2Order : public Interpolator2D
{

public:
    Interpolator2D2Order(Params&, Patch*);
    ~Interpolator2D2Order() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;

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

    void interpolate_em_fields_and_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void interpolate_time_centered_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, int ipart, double* Env_A_abs_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc) override final;

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
