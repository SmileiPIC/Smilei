#ifndef PROJECTORRZ2ORDER_H
#define PROJECTORRZ2ORDER_H

#include <complex>

#include "ProjectorRZ.h"


class ProjectorRZ2Order : public ProjectorRZ {
public:
    ProjectorRZ2Order(Params&, Patch* patch);
    ~ProjectorRZ2Order();

    //! Project global current densities for m=0 (EMfields->Jx_/Jy_/Jz_)
    inline void operator() (std::complex<double>* Jl, std::complex<double>* Jr, std::complex<double>* Jt, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold);

    inline void operator() (std::complex<double>* Jl, std::complex<double>* Jr, std::complex<double>* Jt, Particles &particles, unsigned int ipart,double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold,std::complex<double>* exp_m_theta_old, int imode);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void operator() (std::complex<double>* Jl, std::complex<double>* Jr, std::complex<double>* Jt, std::complex<double>* rho, Particles &particles, unsigned int ipart, double invgf, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold);
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void operator() (std::complex<double>* Jl, std::complex<double>* Jr, std::complex<double>* Jt, std::complex<double>* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim, int* iold, double* deltaold, std::complex<double>* exp_m_theta_old,  int imode);

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void operator() (double* rho, Particles &particles, unsigned int ipart, unsigned int bin, std::vector<unsigned int> &b_dim) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void operator() (Field* Jl, Field* Jr, Field* Jt, Particles &particles, int ipart, LocalFields Jion) override final;

    //!Wrapper
    void operator() (ElectroMagn* EMfields, Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ibin, int clrw, bool diag_flag, bool is_spectral, std::vector<unsigned int> &b_dim, int ispec) override final;

private:
    double one_third;
    // Number of theta modes
    int nmodes;
    double dr;
    double dt;
};

#endif

