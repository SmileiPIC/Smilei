#ifndef PROJECTOR_H
#define PROJECTOR_H

#include "Params.h"
#include "Field.h"

class PicParams;
class Patch;

class ElectroMagn;
class Field;
class Particles;


//----------------------------------------------------------------------------------------------------------------------
//! class Projector: contains the virtual operators used during the current projection
//----------------------------------------------------------------------------------------------------------------------
class Projector {

public:
    //! Creator for the Projector
    Projector(Params&, Patch*);
    virtual ~Projector() {};
    virtual void mv_win(unsigned int shift) = 0;
    virtual void setMvWinLimits(unsigned int shift) = 0;

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    virtual void operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim) = 0;

    //! Project global current densities if Ionization in Species::dynamics,
    virtual void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) = 0;

   //!Wrapper
    virtual void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec) = 0;
private:

};

#endif

