#ifndef PROJECTOR2D_H
#define PROJECTOR2D_H

#include "Projector.h"
#include "PicParams.h"



//----------------------------------------------------------------------------------------------------------------------
//! class Projector2D: defines a virtual method for projection in 1d3v simulations
//----------------------------------------------------------------------------------------------------------------------
class Projector2D : public Projector
{

public:
    //! Constructor for Projector2D
    Projector2D(PicParams& params, SmileiMPI* smpi) : Projector(params, smpi) {};
    virtual ~Projector2D() {};

    virtual void mv_win(unsigned int shift) { i_domain_begin+=shift; }
    virtual void setMvWinLimits(unsigned int shift) {i_domain_begin = shift;}

protected:
    //! Inverse of the spatial step 1/dx
    double dx_inv_;
    double dy_inv_;
    double dx_ov_dt;
    double dy_ov_dt;
    int i_domain_begin;
    int j_domain_begin;
};

#endif

