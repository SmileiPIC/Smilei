#ifndef SPECIESNORM_H
#define SPECIESNORM_H

#include "SpeciesV.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class Params;

//! class SpeciesNorm (Species for which the dynamics is governed by the Lorentz force (Boris pusher))
class SpeciesNorm : public SpeciesV
{

public:
    //! Creator for SpeciesNorm
    SpeciesNorm(Params&, Patch*);
    //! Destructor for SpeciesNorm
    ~SpeciesNorm();

private:

};

#endif
