#ifndef PROJECTOR3D2ORDERSUSCEPTIBILITY_H
#define PROJECTOR3D2ORDERSUSCEPTIBILITY_H

#include "Projector3D.h"


class Projector3D2Order_susceptibility : public Projector3D {
public:
    Projector3D2Order_susceptibility(Params&, Patch* patch);
    ~Projector3D2Order_susceptibility();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    inline void operator() (double* Jx, double* Jy, double* Jz, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void operator() (double* Jx, double* Jy, double* Jz, double* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold);

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void operator() (Field* Jx, Field* Jy, Field* Jz, Particles &particles, int ipart, LocalFields Jion) override final;

    //!Wrapper
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec, int ipart_ref) override final;

    // projects susceptibility  
    void project_susceptibility(double* Chi_envelope, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim, SmileiMPI* smpi, int ithread, double species_mass);

    double dt, dts2, dts4;
    

private:
    double one_third;
};

#endif

