#ifndef PROJECTOR2D4ORDER_H
#define PROJECTOR2D4ORDER_H

#include "Projector2D.h"


class Projector2D4Order : public Projector2D {
public:
    Projector2D4Order(Params&, Patch* patch);
    ~Projector2D4Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    inline void currents(double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void currentsAndDensity(double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, int* iold, double* deltaold);

    //! Project global current charge (EMfields->rho_ , J), for initialization and diags
    void densityFrozen(double* rhoj, Particles &particles, unsigned int ipart, unsigned int type) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents(Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) override final;

    //!Wrapper
    void currentsAndDensityWrapper(ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref = 0) override final;

 private:
    static constexpr double dble_1_ov_384   = 1.0/384.0;
    static constexpr double dble_1_ov_48    = 1.0/48.0;
    static constexpr double dble_1_ov_16    = 1.0/16.0;
    static constexpr double dble_1_ov_12    = 1.0/12.0;
    static constexpr double dble_1_ov_24    = 1.0/24.0;
    static constexpr double dble_19_ov_96   = 19.0/96.0;
    static constexpr double dble_11_ov_24   = 11.0/24.0;
    static constexpr double dble_1_ov_4     = 1.0/4.0;
    static constexpr double dble_1_ov_6     = 1.0/6.0;
    static constexpr double dble_115_ov_192 = 115.0/192.0;
    static constexpr double dble_5_ov_8     = 5.0/8.0;
};

#endif

