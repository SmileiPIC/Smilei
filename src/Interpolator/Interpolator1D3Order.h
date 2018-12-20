#ifndef INTERPOLATOR1D3ORDER_H
#define INTERPOLATOR1D3ORDER_H

#include "Interpolator1D.h"
#include "Field1D.h"

class Interpolator1D3Order : public Interpolator1D {
public:
    Interpolator1D3Order(Params&, Patch*);
    ~Interpolator1D3Order() override final{};
    
    inline void fields    (ElectroMagn* EMfields, Particles &particles, int ipart, int nparts, double* ELoc, double* BLoc);
    void fieldsAndCurrents(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, LocalFields* JLoc, double* RhoLoc) override final;
    void fieldsWrapper     (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int *istart, int *iend, int ithread, int ipart_ref = 0)override final;
    void fieldsSelection (ElectroMagn* EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> * selection) override final;
    void oneField         (Field* field, Particles &particles, int *istart, int *iend, double* FieldLoc) override final;
    
    inline double compute( double* coeff, Field1D* f, int idx) {
        double interp_res =  coeff[0] * (*f)(idx-1)   + coeff[1] * (*f)(idx)   + coeff[2] * (*f)(idx+1) + coeff[3] * (*f)(idx+2);
        return interp_res;
    };
    
protected:
    double dble_1ov6;
    double dble_2ov3;
private:
    inline void coeffs( double xjn ) {
        double xi2, xi3;
        
        // Dual
        id_   = (int)(xjn+0.50);        // position of the 2nd node
        xi  = xjn - (double)id_ + 0.5;    // normalized distance to the central node
        xi2 = xi*xi;
        xi3 = xi*xi2;
        
        // 3rd order interpolation on 4 nodes
        coeffd_[0] = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
        coeffd_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
        coeffd_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
        coeffd_[3]  = xi3*dble_1ov6;
        
        id_ -= index_domain_begin;
        
        // Primal
        ip_ = (int)xjn;            // index of the 2nd node
        xi  = xjn - (double)ip_;   //normalized distance to the 2nd node
        xi2 = xi*xi;
        xi3 = xi2*xi;
        
        // 3rd order interpolation on 4 nodes
        coeffp_[0]  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
        coeffp_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
        coeffp_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
        coeffp_[3]  = xi3*dble_1ov6;
        
        ip_ -= index_domain_begin;
    }
    
    // Last prim index computed
    int ip_;
    // Last dual index computed
    int id_;
    // Last delta computed
    double xi;
    // Interpolation coefficient on Prim grid
    double coeffp_[4];
    // Interpolation coefficient on Dual grid
    double coeffd_[4];

};

#endif
