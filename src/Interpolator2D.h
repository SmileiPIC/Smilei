#ifndef INTERPOLATOR2D_H
#define INTERPOLATOR2D_H


#include "Interpolator.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator 2D
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator2D : public Interpolator
{
public:
	Interpolator2D(PicParams* params, SmileiMPI* smpi): Interpolator(params, smpi) {};
	virtual ~Interpolator2D() {};

	virtual void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc) = 0;
    
protected:
    //! Inverse of the spatial-step
	double dx_inv_;
	double dy_inv_;
	int index_domain_begin;
};

#endif
