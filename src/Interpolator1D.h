#ifndef INTERPOLATOR1D_H
#define INTERPOLATOR1D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 1D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator1D : public Interpolator
{
public:
	Interpolator1D(PicParams* params, SmileiMPI* smpi): Interpolator(params, smpi){;};
	virtual ~Interpolator1D() {};

	virtual void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc) = 0;
    
protected:
    //! Inverse of the spatial-step
	double dx_inv_;
	int index_domain_begin;
};

#endif
