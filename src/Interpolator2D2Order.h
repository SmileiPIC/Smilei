#ifndef INTERPOLATOR2D2ORDER_H
#define INTERPOLATOR2D2ORDER_H


#include "Interpolator2D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D2Order : public Interpolator2D
{

public:
    Interpolator2D2Order(PicParams*, SmileiMPI*);
    ~Interpolator2D2Order();

    void operator() (ElectroMagn* champs, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc);

private:

};//END class

#endif
