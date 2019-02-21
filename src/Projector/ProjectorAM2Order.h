#ifndef PROJECTORAM2ORDER_H
#define PROJECTORAM2ORDER_H

#include <complex>

#include "ProjectorAM.h"


class ProjectorAM2Order : public ProjectorAM
{
public:
    ProjectorAM2Order( Params &, Patch *patch );
    ~ProjectorAM2Order();
    
    //! Project global current densities for m=0 (EMfields->Jl_/Jr_/Jt_)
    inline void currents_mode0( std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold );
    inline void currents( std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, std::complex<double> *exp_m_theta_old, int imode );
    //! Project global current densities (EMfields->Jl_/Jr_/Jt_/rho), diagFields timestep
    inline void currentsAndDensity_mode0( std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, std::complex<double> *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold );
    //! Project global current densities (EMfields->Jl_/Jr_/Jt_/rho), diagFields timestep
    inline void currentsAndDensity( std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, std::complex<double> *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, std::complex<double> *exp_m_theta_old,  int imode );
    
    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplex( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode ) override final;
    
    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final;
    
private:
};

#endif

