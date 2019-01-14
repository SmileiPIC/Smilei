#ifndef PROJECTOR3D_H
#define PROJECTOR3D_H

#include "Projector.h"
#include "Params.h"



//----------------------------------------------------------------------------------------------------------------------
//! class Projector3D: defines a virtual method for projection in 1Dcartesian simulations
//----------------------------------------------------------------------------------------------------------------------
class Projector3D : public Projector
{

public:
    //! Constructor for Projector3D
    Projector3D(Params& params, Patch* patch);
    virtual ~Projector3D() {};

    virtual void mv_win(unsigned int shift) { i_domain_begin+=shift; }
    virtual void setMvWinLimits(unsigned int shift) {i_domain_begin = shift;}

protected:
    //! Inverse of the spatial step 1/dx
    double dx_inv_;
    double dy_inv_;
    double dz_inv_;
    double dx_ov_dt;
    double dy_ov_dt;
    double dz_ov_dt;
    int i_domain_begin;
    int j_domain_begin;
    int k_domain_begin;
    int nscelly;
    int nscellz;
    int nprimy;
    int nprimz;
    int oversize[3];
    double dq_inv[3];
    double *Jx_, *Jy_, *Jz_, *rho_;
    static constexpr double one_third = 1./3.;
};

#endif

