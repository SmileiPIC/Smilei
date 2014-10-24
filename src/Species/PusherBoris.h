/*! @file Pusher.h

 @brief Pusher.h  generic class for the particle pusher

 @author tommaso vinci
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
    PusherBoris(PicParams& params, int ispec);
    ~PusherBoris();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf);

};

#endif

