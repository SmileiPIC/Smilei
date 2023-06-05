/*! @file PusherPonderomotiveBorisBTIS3.h

 @brief PusherPonderomotiveBorisBTIS3.h  generic class for the particle pusher of Boris adapted for envelope model. Only pushes momentum, not their position


 */

#ifndef PUSHERPONDEROMOTIVEBORISBTIS3_H
#define PUSHERPONDEROMOTIVEBORISBTIS3_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPonderomotiveBoris, only pushes momentum of particles interacting with envelope, not their position
//  --------------------------------------------------------------------------------------------------------------------
class PusherPonderomotiveBorisBTIS3 : public Pusher
{
public:
    //! Creator for Pusher
    PusherPonderomotiveBorisBTIS3( Params &params, Species *species );
    ~PusherPonderomotiveBorisBTIS3();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
    
};

#endif
