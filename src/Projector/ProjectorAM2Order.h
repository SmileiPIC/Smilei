#ifndef PROJECTORAM2ORDER_H
#define PROJECTORAM2ORDER_H

#include <complex>

#include "ProjectorAM.h"
#include "ElectroMagnAM.h"


class ProjectorAM2Order : public ProjectorAM
{
public:
    ProjectorAM2Order( Params &, Patch *patch );
    ~ProjectorAM2Order();
    
    inline void currents( ElectroMagnAM *emAM, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, double *array_theta_old, bool diag_flag, int ispec);

    inline void currentsForTasks( ElectroMagnAM *emAM, std::complex<double> *b_Jl, std::complex<double> *b_Jr, std::complex<double> *b_Jt, std::complex<double> *b_rhoAM, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, double *array_theta_old, int bin_shift, int bdim0, bool diag_flag );

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplex( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode ) override final;

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplexOnBuffer( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode, int bdim0, int bin_shift ) override final;

    //! Apply boundary conditions on Rho and J
    void axisBC( ElectroMagnAM *emAM, bool diag_flag ) override final;
    void apply_axisBC(std::complex<double> *rhoj,std::complex<double> *Jl, std::complex<double> *Jr, std::complex<double> *Jt, unsigned int imode, bool diag_flag );
    
    //! Apply boundary conditions on Env_Chi
    void axisBCEnvChi( double *EnvChi ) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    
    //! Project global current densities if Ionization in Species::dynamicsTasks
    void ionizationCurrentsForTasks( double *b_Jx, double *b_Jy, double *b_Jz, Particles &particles, int ipart, LocalFields Jion, int bin_shift ) override final;

    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final;

    //!Wrapper for projection on buffers
    void currentsAndDensityWrapperOnBuffers( double *b_Jx, double *b_Jy, double *b_Jz, double *b_rho, int bin_width, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final {};

    //!Wrapper for projection on AM buffers
    void currentsAndDensityWrapperOnAMBuffers( ElectroMagn *EMfields, std::complex<double> *b_Jl, std::complex<double> *b_Jr, std::complex<double> *b_Jt, std::complex<double> *b_rhoAM, int bin_shift, int bdim0, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, int ipart_ref = 0 ) override final;

    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell = 0, int ipart_ref = 0 ) override final;

    void susceptibilityOnBuffer( ElectroMagn *EMfields, double *b_ChiAM, int bin_shift, int bdim0, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell = 0, int ipart_ref = 0 ) override final;

    

private:
    double dt, dts2, dts4;
};

#endif

