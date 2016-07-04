#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "Params.h"

class Params;
class Patch;
class ElectroMagn;
class Particles;


//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator
{
public:
    Interpolator(Params& params, Patch* patch);
    virtual ~Interpolator() {};
    
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc) = 0;
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread) = 0;
    
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc) = 0;

private:

};//END class

#endif
