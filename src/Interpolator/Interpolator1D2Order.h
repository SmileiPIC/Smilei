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
    void operator() (Field* field, Particles &particles, int *istart, int *iend, double* FieldLoc) override final;
    
    inline double compute( double* coeff, Field1D* f, int idx) {
        double interp_res =  coeff[0] * (*f)(idx-1)   + coeff[1] * (*f)(idx)   + coeff[2] * (*f)(idx+1);
        return interp_res;
    };
    
private:
    inline void coeffs( double xjn ) {
        double xjmxi2;
        
        // Dual
        id_      = round(xjn+0.5);        // index of the central point
        xjmxi  = xjn - (double)id_ +0.5;  // normalized distance to the central node
        xjmxi2 = xjmxi*xjmxi;            // square of the normalized distance to the central node
        
        // 2nd order interpolation on 3 nodes
        coeffd_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
        coeffd_[1] = (0.75-xjmxi2);
        coeffd_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
        
        id_ -= index_domain_begin;
        
        // Primal
        ip_      = round(xjn);      // index of the central point
        xjmxi  = xjn -(double)ip_;  // normalized distance to the central node
        xjmxi2 = pow(xjmxi,2);      // square of the normalized distance to the central node
        
        // 2nd order interpolation on 3 nodes
        coeffp_[0] = 0.5 * (xjmxi2-xjmxi+0.25);
        coeffp_[1] = (0.75-xjmxi2);
        coeffp_[2] = 0.5 * (xjmxi2+xjmxi+0.25);
        
        ip_ -= index_domain_begin;
    }
    
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xjmxi;
    // Interpolation coefficient on Prim grid
    double coeffp_[3];
    // Interpolation coefficient on Dual grid
    double coeffd_[3];
    
    
};//END class

#endif
