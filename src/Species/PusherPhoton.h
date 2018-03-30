/* --------------------------------------------------------------------------------------------------------------------

! @file PusherPhoton.h

 @brief PusherPhoton.h  generic class for the photon pusher.

 @date 2017-07-21

 --------------------------------------------------------------------------------------------------------------------
 */

#ifndef PUSHERPHOTON_H
#define PUSHERPHOTON_H

#include "Pusher.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherPhoton
//  --------------------------------------------------------------------------------------------------------------------
class PusherPhoton : public Pusher {
public:
    //! Creator for Pusher
    PusherPhoton(Params& params, Species *species);
    ~PusherPhoton();
    //! Overloading of () operator
    virtual void operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ipart_ref = 0);

};

#endif
