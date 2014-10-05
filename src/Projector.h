#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "PicParams.h"
#include "SmileiMPI.h"
#include "Field.h"

class ElectroMagn;
class Field;
class Particles;


//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector {

public:
    //! Creator for the Projector
    Projector(PicParams&, SmileiMPI*) {};
    virtual ~Projector() {};
    virtual void mv_win(unsigned int shift) = 0;
    virtual void setMvWinLimits(unsigned int shift) = 0;

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf) = 0;

    // Projection by species, in Species::dynamics, ElectroMagn::initRhoJ
    virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf) = 0;

    //! Project global current charge (EMfields->rho_)
    //! Used in Species::dynamics if time_frozen
    virtual void operator() (Field* rho, Particles &particles, int ipart) = 0;


    //! Project local current densities if particles sorting activated in Species::dynamics
    virtual void operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0) = 0;

    //! Project global current densities if Ionization in Species::dynamics,
    virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) = 0;

private:

};

#endif

