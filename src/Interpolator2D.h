#ifndef INTERPOLATOR2D_H
#define INTERPOLATOR2D_H


#include "Interpolator.h"
#include "SmileiMPI_Cart2D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 2D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D : public Interpolator
{
public:
    Interpolator2D(PicParams&params, SmileiMPI*smpi): Interpolator(params, smpi) {
        SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
        i_domain_begin = smpi2D->getCellStartingGlobalIndex(0);
        j_domain_begin = smpi2D->getCellStartingGlobalIndex(1);
    };

    virtual ~Interpolator2D() {};
    virtual void mv_win(unsigned int shift) {i_domain_begin += shift;}
    virtual void setMvWinLimits(unsigned int shift) {i_domain_begin = shift;}

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

protected:
    //! Inverse of the spatial-step
    double dx_inv_;
    double dy_inv_;
    int i_domain_begin;
    int j_domain_begin;
};

#endif
