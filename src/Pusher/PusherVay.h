/*! @file PusherVay.h

 @brief PusherVay.h  generic class for the particle pusher of J.L. Vay.

 #details The description of the J.L. Vay pusher can be found in this reference:
          http://dx.doi.org/10.1063/1.2837054

 @date 2017-04-10
 */

#ifndef PUSHERVAY_H
#define PUSHERVAY_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherVay
//  --------------------------------------------------------------------------------------------------------------------
class PusherVay : public Pusher
{
public:
    //! Creator for Pusher
    PusherVay( Params &params, Species *species );
    ~PusherVay();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
    
};

#endif
