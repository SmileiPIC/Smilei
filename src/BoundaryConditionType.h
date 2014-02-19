/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include "Particles.h"

inline void refl_particle( Particles &particles, int ipart, int direction, double limit_pos ) {
    particles.position(direction, ipart) = limit_pos - particles.position(direction, ipart);
    particles.momentum(direction, ipart) = -particles.momentum(direction, ipart);

}

#endif

