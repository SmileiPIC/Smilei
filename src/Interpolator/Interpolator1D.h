#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 1D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D : public Interpolator
{
public:
    Interpolator1D(Params &params, Patch *patch) ;

    virtual ~Interpolator1D() {};
    virtual void mv_win(unsigned int shift) {index_domain_begin += shift;}
    virtual void setMvWinLimits(unsigned int shift) {index_domain_begin = shift;}

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    unsigned int index_domain_begin;
};

#endif
