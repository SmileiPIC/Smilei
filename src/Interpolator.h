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
    Interpolator(PicParams*, SmileiMPI*) {};
    virtual ~Interpolator() {};

    virtual void operator() (ElectroMagn* champs, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;

private:

};//END class

#endif
