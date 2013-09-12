/*! @file Pusher.h
 
  @brief Pusher.h  generic class for the particle pusher
 
  @author tommaso vinci
  @date 2013-02-15
*/

#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "PicParams.h"

#include "Particle.h"

class SmileiMPI;

class PartBoundCond {
public:
    PartBoundCond( PicParams *params, int ispec, SmileiMPI* smpi );
    ~PartBoundCond();

    void (*bc_west)  ( Particle* part, int direction, double limit_pos );
    void (*bc_east)  ( Particle* part, int direction, double limit_pos );
    void (*bc_south) ( Particle* part, int direction, double limit_pos );
    void (*bc_north) ( Particle* part, int direction, double limit_pos );
    void (*bc_bottom)( Particle* part, int direction, double limit_pos );
    void (*bc_up)    ( Particle* part, int direction, double limit_pos );

    inline int apply( Particle* part ) {

    	if ( part->position(0) <  x_min ) {
    		if (bc_west==NULL) return 0;
    		else { (*bc_west)( part, 0, 2.*x_min ); return 1; }
    	}
    	else if ( part->position(0) >= x_max ) {
    		if (bc_east==NULL) return 0;
    		else { (*bc_east)( part, 0, 2.*x_max ); return 1; }
    	}
        else if (nDim_particle == 2) {

			if ( part->position(1) <  y_min ) {
				if (bc_south==NULL) return 0;
				else { (*bc_south)( part, 1, 2.*y_min ); return 1; }
			}
			else if ( part->position(1) >= y_max ) {
				if (bc_north==NULL) return 0;
				else { (*bc_north)( part, 1, 2.*y_max ); return 1; }
			}

			else if (nDim_particle == 3) {

				if ( part->position(2) <  z_min ) {
					if (bc_bottom==NULL) return 0;
					else { (*bc_bottom)( part, 2, 2.*z_min ); return 1; }
				}
				else if ( part->position(2) >= z_max ) {
					if (bc_up==NULL) return 0;
					else { (*bc_up)( part, 2, 2.*z_max ); return 1; }
				}
			} // end if (nDim_particle == 2)
        } // end if (nDim_particle == 3)

    	return 1;
    };

 private:
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    double z_min;
    double z_max;

    int nDim_particle;

};

#endif

