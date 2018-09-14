#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"

#include "Field1D.h"
//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1Dcartesian simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order : public Interpolator1D
{

public:
    Interpolator1D2Order(Params&, Patch*);
    ~Interpolator1D2Order() override final{};
    
    inline void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0) override final;
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final;
    void operator() (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;
    
    inline double compute( double* coeff, Field1D* f, int idx) {
        double interp_res =  coeff[0] * (*f)(idx-1)   + coeff[1] * (*f)(idx)   + coeff[2] * (*f)(idx+1);
        return interp_res;
    };
    
    void interpolate_em_fields_and_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void interpolate_envelope_and_old_envelope( ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0 ) override final;
    void interpolate_envelope_and_susceptibility(ElectroMagn* EMfields, Particles &particles, int ipart, double* Env_A_abs_Loc, double* Env_Chi_Loc, double* Env_E_abs_Loc) override final;

private:
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xjmxi,deltax;
    // Interpolation coefficient on Prim grid
    double coeffp_[3];
    // Interpolation coefficient on Dual grid
    double coeffd_[3];
    
    
};//END class

#endif
