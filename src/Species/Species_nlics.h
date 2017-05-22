// ---------------------------------------------------------------------------------------------------------------------
//
//! \file Species_nlics.h
//
//! \brief Species_nlics.h  generic class for the species that use the Monte-Carlo
//!  pusher for the Non Linear Inverse Compton Scattering.
//
//! \date 2017-05-12
//
// ---------------------------------------------------------------------------------------------------------------------

#ifndef SPECIESNLICS_H
#define SPECIESNLICS_H

#include "Species.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class Params;

//! class Species_nlics
//! Species for which the dynamics is governed by the Lorentz force and the
//! Nonlinear Inverse Compton Scattering quantum radiation reaction force
//! (Monte-Carlo process)
class Species_nlics : public Species
{
public:
    //! Creator for Species_nlics
    Species_nlics(Params&, Patch*);
    //! Destructor for Species_nlics
    ~Species_nlics();

private:

};

#endif
