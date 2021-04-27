/*! @file PusherBorisNR.h

 @brief PusherBorisNR.h  Non-Relativistic Boris Particle Pusher

 @date 2016-04-11
 */

#ifndef PUSHERBORISNR_H
#define PUSHERBORISNR_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBorisNR
//  --------------------------------------------------------------------------------------------------------------------
class PusherBorisNR : public Pusher
{
public:
    //! Creator for Pusher
    PusherBorisNR( Params &params, Species *species );
    ~PusherBorisNR();
    
    //! Overriding operator()
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
    
};

#endif
