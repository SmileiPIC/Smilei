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
    
    //! Overloading of () operator
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Bpart, double& gf) = 0;

    //\todo Why not put this in the Species class? (MG to JD)
    
    //! Method used to get the particle mass
	inline double getMass()   {return mass_  ;};
    
    //
	inline double getCharge() {return charge_;};

protected:
	double dt, dts2;
	double mass_;
	double charge_;
	double charge_over_mass_;
    
};//END class

#endif

