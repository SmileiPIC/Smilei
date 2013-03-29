
#ifndef PUSHERSKLV_H
#define PUSHERSKLV_H

#include "PusherSklv.h"
#include <iostream>

class Particle;

class PusherSklv : public Pusher {
public:
	PusherSklv(PicParams *params, int ispec);
	virtual void operator() (Particle* part, chLocaux Epart, chLocaux Bpart, double& gf) { std::cout << "\tSokolov Push particle" << std::endl; };
private:
};

#endif

