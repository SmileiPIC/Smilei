#ifndef PUSHERSKLV_H
#define PUSHERSKLV_H

#include <iostream>

#include "PusherSklv.h"

class Particles;

//  --------------------------------------------------------------------------------------------------------------------
//! Class PusherSklv (not implemented for now)
//  --------------------------------------------------------------------------------------------------------------------
class PusherSklv : public Pusher {
public:
    //! Creator for PusherSklv
    PusherSklv(PicParams& params, int ispec);
    //! Overloading of () operator
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf) {
        std::cout << "\tSokolov Push particle" << std::endl;
    };
private:

};

#endif

