#ifndef SPECIESNORM_H
#define SPECIESNORM_H

#include "Species.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PicParams;

//! class Species_norm (Species for which the dynamics is governed by the Lorentz force (Boris pusher))
class Species_norm : public Species
{

public:
    //! Creator for Species_norm
    Species_norm(PicParams&, int, SmileiMPI*);
    //! Destructor for Species_norm
    ~Species_norm();

private:

};

#endif
