#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order : public Interpolator1D
{

public:
    Interpolator1D2Order(PicParams*, SmileiMPI*);
    ~Interpolator1D2Order();

    void operator() (ElectroMagn* champs, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);

private:

};//END class

#endif
