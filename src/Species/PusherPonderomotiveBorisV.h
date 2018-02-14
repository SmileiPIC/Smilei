/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher of Boris adapted for envelope model. Only pushes momentum, not their position

 @date 2013-02-15
 */

#ifndef PUSHERPONDEROMOTIVEBORISV_H
#define PUSHERPONDEROMOTIVEBORISV_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPonderomotiveBorisV, only pushes momentum of particles interacting with envelope, not their position
//  --------------------------------------------------------------------------------------------------------------------
class PusherPonderomotiveBorisV : public Pusher {
public:
    //! Creator for Pusher
    PusherPonderomotiveBorisV(Params& params, Species *species);
    ~PusherPonderomotiveBorisV();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread);

};

#endif

