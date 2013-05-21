#ifndef INTERPOLATOR1D2ORDER_H
#define INTERPOLATOR1D2ORDER_H


#include "Interpolator1D.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Class for 2nd order interpolator for 1d3v simulations
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D2Order : public Interpolator1D
{
    
public:
    //! Creator for Interpolator1D2Order
	Interpolator1D2Order(PicParams *);
    
    // \todo Comment
    //! Overloading of the () operator
	void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc);
    
private:
    
};//END class

#endif
