#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
    Projector1D2Order(PicParams&, SmileiMPI* smpi);
    ~Projector1D2Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! Not used for now
    void operator() (ElectroMagn* EMfields, Particles &particles, int ipart, double gf);

    // Projection by species, in Species::dynamics, ElectroMagn::initRhoJ
    void operator() (Field* Jx, Field* Jy, Field* Jz, Field* rho, Particles &particles, int ipart, double gf);

    //! Project global current charge (EMfields->rho_)
    //! Used in Species::dynamics if time_frozen
    void operator() (Field* rho, Particles &particles, int ipart);

     //! Project local current densities if particles sorting activated in Species::dynamics
    void operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, int ipart, double gf, unsigned int bin, unsigned int b_dim0);

    //! Project global current densities if Ionization in Species::dynamics,
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion);

private:
    double dx_ov_dt;
};

#endif

