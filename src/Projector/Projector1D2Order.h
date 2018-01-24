#ifndef PROJECTOR1D2ORDER_H
#define PROJECTOR1D2ORDER_H

#include "Projector1D.h"

class Projector1D2Order : public Projector1D {
public:
    Projector1D2Order(Params&, Patch* patch);
    ~Projector1D2Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    inline void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* delta);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* delta);

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) override final;

    //!Wrapper
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec) override final;

private:
    double dx_ov_dt;
};

#endif

