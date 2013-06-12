/*! @file Pusher.h
 
 @brief Pusher.h  generic class for the particle pusher
 
 @author tommaso vinci
 @date 2013-02-15
 */

#ifndef PUSHERBORIS_H
#define PUSHERBORIS_H

#include "Pusher.h"

class PusherBoris : public Pusher {
public:
    PusherBoris(PicParams *params, int ispec);
    ~PusherBoris();
	virtual void operator() (Particle* part, LocalFields Epart, LocalFields Bpart, double& gf);
    
};

#endif

