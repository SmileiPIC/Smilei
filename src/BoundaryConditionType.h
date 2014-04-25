/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include "Particles.h"

//!
//! int function( Particles &particles, int ipart, int direction, double limit_pos )
//!     returns :
//!         0 if particle ipart have ti be deleted from current process (MPI or BC)
//!         1 otherwise
//!

inline int refl_particle( Particles &particles, int ipart, int direction, double limit_pos ) {
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);
    return 1;

}

inline int supp_particle( Particles &particles, int ipart, int direction, double limit_pos ) {
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    particles.charge(ipart) = 0;
    return 0;
}

inline int stop_particle( Particles &particles, int ipart, int direction, double limit_pos ) {
    particles.position(direction, ipart) = particles.position_old(direction, ipart);
    return 0;

}


#endif
