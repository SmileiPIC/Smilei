#ifndef SPECIESNORMV_H
#define SPECIESNORMV_H

#include "SpeciesV.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class Params;

//! class SpeciesNorm (Species for which the dynamics is governed by the Lorentz force (Boris pusher))
class SpeciesNormV : public SpeciesV
{

public:
    //! Creator for SpeciesNorm
    SpeciesNormV( Params &, Patch * );
    //! Destructor for SpeciesNorm
    ~SpeciesNormV();
    
private:

};

#endif
