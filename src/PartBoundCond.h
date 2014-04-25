/*! @file Pusher.h

  @brief Pusher.h  generic class for the particle pusher

  @author tommaso vinci
  @date 2013-02-15
*/

#ifndef PARTBOUNDCOND_H
#define PARTBOUNDCOND_H

#include "PicParams.h"
#include "Particles.h"

class SmileiMPI;

class PartBoundCond {
public:
    PartBoundCond( PicParams *params, int ispec, SmileiMPI* smpi );
    ~PartBoundCond();

    int (*bc_west)  ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_east)  ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_south) ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_north) ( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_bottom)( Particles &particles, int ipart, int direction, double limit_pos );
    int (*bc_up)    ( Particles &particles, int ipart, int direction, double limit_pos );

    inline int apply( Particles &particles, int ipart ) {

        if ( particles.position(0, ipart) <  x_min ) {
            if (bc_west==NULL) return 0;
            else {
                return (*bc_west)( particles, ipart, 0, 2.*x_min );
            }
        }
        else if ( particles.position(0, ipart) >= x_max ) {
            if (bc_east==NULL) return 0;
            else {
                return (*bc_east)( particles, ipart, 0, 2.*x_max );
            }
        }
        else if (nDim_particle == 2) {

            if ( particles.position(1, ipart) <  y_min ) {
                if (bc_south==NULL) return 0;
                else {
                    return (*bc_south)( particles, ipart, 1, 2.*y_min );
                }
            }
            else if ( particles.position(1, ipart) >= y_max ) {
                if (bc_north==NULL) return 0;
                else {
                    return (*bc_north)( particles, ipart, 1, 2.*y_max );
                }
            }

            else if (nDim_particle == 3) {

                if ( particles.position(2, ipart) <  z_min ) {
                    if (bc_bottom==NULL) return 0;
                    else {
                        return (*bc_bottom)( particles, ipart, 2, 2.*z_min );
                    }
                }
                else if ( particles.position(2, ipart) >= z_max ) {
                    if (bc_up==NULL) return 0;
                    else {
                        return (*bc_up)( particles, ipart, 2, 2.*z_max );
                    }
                }
            } // end if (nDim_particle == 2)
        } // end if (nDim_particle == 3)

        return 1;
    };

    void moveWindow_x(double shift);

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

