#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "PicParams.h"
#include "SmileiMPI.h"


class ElectroMagn;
class Particles;



//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator
{
public:
    Interpolator(PicParams&, SmileiMPI*) {};
    virtual ~Interpolator() {};
    virtual void mv_win(unsigned int shift) = 0;
    virtual void setMvWinLimits(unsigned int shift) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;

    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

private:

};//END class

#endif
