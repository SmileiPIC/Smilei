/*! @file Pusher.h
 
 @brief Pusher.h  generic class for the particle pusher
 
 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef BOUNDARYCONDITIONTYPE_H
#define BOUNDARYCONDITIONTYPE_H

#include "Particle.h"

inline void refl_particle( Particle* part, double limit_pos ) {
        part->position(0) = limit_pos - part->position(0);
        part->momentum(0) = -part->momentum(0);

}

#endif

