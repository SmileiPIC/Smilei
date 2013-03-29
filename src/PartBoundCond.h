/*! @file Pusher.h
 
 @brief Pusher.h  generic class for the particle pusher
 
 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "PicParams.h"
#include "Particle.h"

class PartBoundCond {
public:
	PartBoundCond( PicParams *params, int ispec);
	~PartBoundCond();

	inline int locateParticle( Particle* part ) {               
        	int place = 0;

		//! \todo{locate : define 0 to ...}
	        if      ( part->position(0) <= x_min) return 1;
	       	else if ( part->position(0) >= x_max) return 2;
		else if (bc_north  != NULL) {
		        if ( part->position(1) <= y_min)      return 3;
	        	else if ( part->position(1) >= y_max) return 4;
	        	else if ( part->position(2) <= z_min) return 5;
		        else if ( part->position(2) >= z_max) return 6;
		}
	        return  place;
	};

	void (*bc_east)  ( Particle* part, double limit_pos );
	void (*bc_west)  ( Particle* part, double limit_pos );
	void (*bc_north) ( Particle* part, double limit_pos );
	void (*bc_south) ( Particle* part, double limit_pos );
	void (*bc_bottom)( Particle* part, double limit_pos );
	void (*bc_up)    ( Particle* part, double limit_pos );

	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;

};

#endif

