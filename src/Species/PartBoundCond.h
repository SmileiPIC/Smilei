/*! @file Pusher.h

  @brief Pusher.h  generic class for the particle pusher

  @author tommaso vinci
  @author Arnaud Beck
  @date 2013-02-15
*/

#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "PicParams.h"
#include "Particles.h"

class SmileiMPI;

class PartBoundCond {
public:
    PartBoundCond( PicParams& params, int ispec, SmileiMPI* smpi );
    ~PartBoundCond();

    int (*bc_west)  ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_east)  ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_south) ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_north) ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_bottom)( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_up)    ( Particles &particles, int ipart, int direction, double limit_pos );

    //Here particles boundary conditions are applied. Conditions along X are applied first, then Y, then Z.
    //The decision whether the particle is added or not on the Exchange Particle List is defined by the final
    //value of keep_part. 
    //Be careful, once an a BC along a given dimension set keep_part to 0, it will remain to 0. 
    //int keep_part; // 0 if particle leave the proc, 1 if particle is kept.
    inline int apply( Particles &particles, int ipart ) {

        int keep_part = 1;
        if ( particles.position(0, ipart) <  x_min ) {
            if (bc_west==NULL) keep_part = 0;
            else {
                keep_part = (*bc_west)( particles, ipart, 0, 2.*x_min );
            }
        }
        else if ( particles.position(0, ipart) >= x_max ) {
            if (bc_east==NULL) keep_part = 0;
            else {
                keep_part = (*bc_east)( particles, ipart, 0, 2.*x_max );
            }
        }
        if (nDim_particle >= 2) {

            if ( particles.position(1, ipart) <  y_min ) {
		if (bc_south==NULL) keep_part = 0;
                else {
                    keep_part *= (*bc_south)( particles, ipart, 1, 2.*y_min );
                }
            }
            else if ( particles.position(1, ipart) >= y_max ) {
		if (bc_north==NULL) keep_part = 0;
                else {
                    keep_part *= (*bc_north)( particles, ipart, 1, 2.*y_max );
                }
            }

            if (nDim_particle == 3) {

                if ( particles.position(2, ipart) <  z_min ) {
                    if (bc_bottom==NULL) keep_part = 0;
                    else {
                        keep_part *= (*bc_bottom)( particles, ipart, 2, 2.*z_min );
                    }
                }
                else if ( particles.position(2, ipart) >= z_max ) {
                    if (bc_up==NULL) keep_part = 0;
                    else {
                        keep_part *= (*bc_up)( particles, ipart, 2, 2.*z_max );
                    }
                }
            } // end if (nDim_particle == 3)
        } // end if (nDim_particle >= 2)

        return keep_part;
    };

    void moveWindow_x(double shift, SmileiMPI* smpi );

    inline void updateMvWinLimits( double x_moved ) {
	x_min += x_moved;
	x_max += x_moved;
    }

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

