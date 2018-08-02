#ifndef INTERPOLATOR3D2ORDER_ENV_H
#define INTERPOLATOR3D2ORDER_ENV_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D2Order_env : public Interpolator3D
{

public:
    Interpolator3D2Order_env(Params&, Patch*);
    ~Interpolator3D2Order_env() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc) {};
    //inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc, double* PHILoc, double* GradPHILoc);
    //inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* PHILoc, double* GradPHILoc, double* PHIoldLoc, double* GradPHIoldLoc, bool prevent_error_overloading);

    void interpolate_em_fields_and_envelope(ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc, double* PHILoc, double* GradPHILoc);
    void interpolate_envelope_and_old_envelope(ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* PHILoc, double* GradPHILoc, double* PHIoldLoc, double* GradPHIoldLoc);
    void interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, int ipart, double* Env_A_abs_Loc, double* Env_Ar_Loc, double* Env_Ai_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc);

    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final ;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;

    inline double compute( double* coeffx, double* coeffy, double* coeffz, Field3D* f, int idx, int idy, int idz) {
	double interp_res(0.);
        //unroll ?
	for (int iloc=-1 ; iloc<2 ; iloc++) {
	    for (int jloc=-1 ; jloc<2 ; jloc++) {
                for (int kloc=-1 ; kloc<2 ; kloc++) {
                    interp_res += *(coeffx+iloc) * *(coeffy+jloc) * *(coeffz+kloc) * (*f)(idx+iloc,idy+jloc,idz+kloc);
                }
	    }
	}
	return interp_res;
    };  

    // Last prim index computed
    int ip_, jp_, kp_;
       // Last dual index computed
    int id_, jd_, kd_;

    // Last delta computed
    double deltax, deltay, deltaz;

private:

    // Interpolation coefficient on Prim grid
    double coeffxp_[3], coeffyp_[3], coeffzp_[3];
    // Interpolation coefficient on Dual grid
    double coeffxd_[3], coeffyd_[3], coeffzd_[3];


};//END class

#endif
