#ifndef PUSHER_H
#define PUSHER_H


#include "PicParams.h"
#include "Field.h"


class Particle;


//  --------------------------------------------------------------------------------------------------------------------
//! Class Pusher
//  --------------------------------------------------------------------------------------------------------------------
class Pusher
{

public:
    //! Creator for Pusher
	Pusher(PicParams *params, int ispec);
	virtual ~Pusher();
    
    //! Overloading of () operator
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Bpart, double& gf) = 0;

    //! Method used to get the particle mass
	//inline double getMass()   {return mass_  ;};
    
    //! Method used to get the particle charge
	//inline double getCharge() {return charge_;};

protected:
	double dt, dts2;
	// mass_ and charge_ relative to Species but used in the particle pusher
	double mass_;
	double one_over_mass_;

	int nDim_;
    
};//END class

#endif

