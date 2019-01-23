/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERBORISV_H
#define PUSHERBORISV_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBorisV
//  --------------------------------------------------------------------------------------------------------------------
class PusherBorisV : public Pusher
{
public:
    //! Creator for Pusher
    PusherBorisV( Params &params, Species *species );
    ~PusherBorisV();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref = 0 );
    
};

#endif

