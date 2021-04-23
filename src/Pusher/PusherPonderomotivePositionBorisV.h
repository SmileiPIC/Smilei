/*! @file PusherPonderomotivePositionBorisV.h

generic class for the particle pusher of Boris adapted for envelope model. Only pushes position, not their momentum


 */

#ifndef PUSHERPONDEROMOTIVEPOSITIONBORISV_H
#define PUSHERPONDEROMOTIVEPOSITIONBORISV_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPonderomotiveBorisV, only pushes momentum of particles interacting with envelope, not their position
//  --------------------------------------------------------------------------------------------------------------------
class PusherPonderomotivePositionBorisV : public Pusher
{
public:
    //! Creator for Pusher
    PusherPonderomotivePositionBorisV( Params &params, Species *species );
    ~PusherPonderomotivePositionBorisV();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
    
};

#endif
