/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @date 2013-02-15
 */

#ifndef PUSHERBORIS_H
#define PUSHERBORIS_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBorisV
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris : public Pusher
{
public:
    //! Creator for Pusher
    PusherBoris( Params &params, Species *species );
    ~PusherBoris();
    //! Overloading of () operator
    void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 ) override;
};

#endif
