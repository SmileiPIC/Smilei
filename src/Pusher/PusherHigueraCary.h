/*! @file PusherHigueraCary.h

  @brief PusherHigueraCary.h  generic class for the particle pusher of A.V. Higuera and J.R. Cari.

  @details See article https://arxiv.org/abs/1701.05605

  @date 2017-03-31
 */

#ifndef PUSHERHIGUERACARY_H
#define PUSHERHIGUERACARY_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherHigueraCary
//  --------------------------------------------------------------------------------------------------------------------
class PusherHigueraCary : public Pusher
{
public:
    //! Creator for Pusher
    PusherHigueraCary( Params &params, Species *species );
    ~PusherHigueraCary();
    //! Overloading of () operator
    virtual void operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset = 0 );
};

#endif
