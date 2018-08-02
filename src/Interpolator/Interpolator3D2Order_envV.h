#ifndef INTERPOLATOR3D2ORDER_ENVV_H
#define INTERPOLATOR3D2ORDER_ENVV_H


#include "Interpolator3D.h"
#include "Field3D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator3D2Order_envV : public Interpolator3D
{

public:
    Interpolator3D2Order_envV(Params&, Patch*);
    ~Interpolator3D2Order_envV() override final {};

    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double* ELoc, double* BLoc);
    // Sorting
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final;
    // Probes
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final {};

    void interpolate_em_fields_and_envelope(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref);
    void interpolate_envelope_and_old_envelope(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref);
    void interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, double* Env_A_abs_Loc, double* Env_Ar_Loc, double* Env_Ai_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc);

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

private:


};//END class

#endif
