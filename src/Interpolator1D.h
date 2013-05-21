#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 1D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D : public Interpolator
{
public:
    //! Creator for Interpolator1D
	Interpolator1D(PicParams * params): Interpolator(params){;};
    
    //! \todo Idem why is () overloaded?
    //! Overloading of the () operator
	virtual void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc) = 0;
    
protected:
    //! Inverse of the spatial-step
	double dx_inv_;
    
};//END class

#endif
