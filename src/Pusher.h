/*! @file Pusher.h
 
 @brief Pusher.h  generic class for the particle pusher
 
 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef PUSHER_H
#define PUSHER_H

#include "PicParams.h"
#include "Field.h"

class Particle;

class Pusher {
public:
	Pusher(PicParams *params, int ispec);
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Bpart, double& gf) = 0;

	inline double getMass()   {return mass_  ;};
	inline double getCharge() {return charge_;};

protected:
	double dt, dts2;
	double mass_;
	double charge_;
	double charge_over_mass_;
    
};

#endif

