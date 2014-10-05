#ifndef PUSHERSKLV_H
#define PUSHERSKLV_H

#include <iostream>

#include "PusherSklv.h"

class Particles;

class PusherSklv : public Pusher {
public:
    PusherSklv(PicParams& params, int ispec);
    virtual void operator() (Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double& gf) {
        std::cout << "\tSokolov Push particle" << std::endl;
    };
private:

};

#endif

