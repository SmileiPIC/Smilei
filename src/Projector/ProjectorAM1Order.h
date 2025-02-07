#ifndef PROJECTORAM1ORDER_H
#define PROJECTORAM1ORDER_H

#include <complex>

#include "ProjectorAM.h"
#include "ElectroMagnAM.h"


class ProjectorAM1Order : public ProjectorAM
{
public:
    ProjectorAM1Order( Params &, Patch *patch );
    ~ProjectorAM1Order();
    
    inline void currents( ElectroMagnAM *emAM, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, bool diag_flag, int ispec);
    
    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplex( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode ) override final;

    //! Apply boundary conditions on Rho and J
    void axisBC( ElectroMagnAM *emAM, bool diag_flag ) override final;
    void apply_axisBC(std::complex<double> *rho, unsigned int imode, unsigned int nonzeromode);

    //! Apply boundary conditions on Env_Chi
    void axisBCEnvChi( double *EnvChi ) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;

    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final;

    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell = 0, int ipart_ref = 0 ) override final;

private:
    double dt, dts2, dts4;
};

#endif
