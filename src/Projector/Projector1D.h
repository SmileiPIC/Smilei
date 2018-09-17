#ifndef PROJECTOR1D_H
#define PROJECTOR1D_H

#include "Projector.h"
#include "Params.h"



//----------------------------------------------------------------------------------------------------------------------
//! class Projector1D: defines a virtual method for projection in 1Dcartesian simulations
//----------------------------------------------------------------------------------------------------------------------
class Projector1D : public Projector
{

public:
    //! Constructor for Projector1D
    Projector1D(Params& params, Patch* patch);
    virtual ~Projector1D() {};
    virtual void mv_win(unsigned int shift) { index_domain_begin+=shift; }
    virtual void setMvWinLimits(unsigned int shift) {index_domain_begin = shift;}

protected:
    //! Inverse of the spatial step 1/dx
    double dx_inv_;
    int index_domain_begin;
    double *Jx_, *Jy_, *Jz_, *rho_;

private:

};

#endif

