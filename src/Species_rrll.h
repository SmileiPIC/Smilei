#ifndef SPECIESRRLL_H
#define SPECIESRRLL_H

#include "Species.h"
#include <string>


class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PicParams;

//! class Species_rrLL
//! Species for which the dynamics is governed by the Lorentz force and the classical radiation reaction force:
//! Boris pusher + first order spliting
class Species_rrLL : public Species
{
    
public:
    //! Creator for Species_rrLL
    Species_rrLL(PicParams*, unsigned int);
    //! Destructor for Species_rrLL
    ~Species_rrLL();
    
    
private:
    
};

#endif
