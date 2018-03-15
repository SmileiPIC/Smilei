/*! @file PusherBoris.h

 @brief PusherBoris.h  generic class for the particle pusher of Boris.

 @date 2013-02-15
 */

#ifndef PUSHERBORIS_H
#define PUSHERBORIS_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherBoris
//  --------------------------------------------------------------------------------------------------------------------
class PusherBoris : public Pusher {
public:
    //! Creator for Pusher
    PusherBoris(Params& params, Species *species);
    ~PusherBoris();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ipart_ref = 0);

};

#endif

