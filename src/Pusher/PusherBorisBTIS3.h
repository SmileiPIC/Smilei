/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERBORISBTIS3_H
#define PUSHERBORISBTIS3_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBorisV
//  --------------------------------------------------------------------------------------------------------------------
class PusherBorisBTIS3 : public Pusher
{
public:
    //! Creator for Pusher
    PusherBorisBTIS3( Params &params, Species *species );
    ~PusherBorisBTIS3();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
    
};

#endif
