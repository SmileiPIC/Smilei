#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include "Field.h"
#include "PicParams.h"
#include "SmileiMPI.h"


class ElectroMagn;
class Particle;



//  --------------------------------------------------------------------------------------------------------------------
//! Class Interpolator
//  --------------------------------------------------------------------------------------------------------------------
class Interpolator
{
public:
    //! Creator for Interpolator
	Interpolator(PicParams *){;};
    
    //! \todo Why is the () operator overloaded? (MG to JD)
    //! Overloading of the () operator
	virtual void operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc) = 0;
private:
    
};//END class

#endif
