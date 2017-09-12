#ifndef PROJECTORRZ_H
#define PROJECTORRZ_H

#include "Projector.h"
#include "Params.h"



//----------------------------------------------------------------------------------------------------------------------
//! class ProjectorRZ: defines a virtual method for projection in 1d3v simulations
//----------------------------------------------------------------------------------------------------------------------
class ProjectorRZ : public Projector
{

public:
    //! Constructor for ProjectorRZ
    ProjectorRZ(Params& params, Patch* patch);
    virtual ~ProjectorRZ() {};

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

