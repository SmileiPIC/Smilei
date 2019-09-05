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
    
    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basicForComplex( std::complex<double> *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int imode ) override final;
    
    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final;

    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell = 0, int ipart_ref = 0 ) override final;
    
private:
};

#endif

