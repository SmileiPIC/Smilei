#ifndef PROJECTORAM2ORDERV_H
#define PROJECTORAM2ORDERV_H

#include "ProjectorAM.h"


class ProjectorAM2OrderV : public ProjectorAM
{
public:
    ProjectorAM2OrderV( Params &, Patch *patch );
    ~ProjectorAM2OrderV();
    
    //! Project global current densities (EMfields->Jl_/Jr_/Jt_)
    void currents(ElectroMagnAM *emAM, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, std::complex<double> *array_eitheta_old, int ipart_ref = 0 );
    ////! Project global current densities (EMfields->Jl_/Jr_/Jt_/rho), diagFields timestep
    //inline void currentsAndDensity( double *Jl, double *Jr, double *Jt, double *rho, Particles &particles, unsigned int istart, unsigned int iend, std::vector<double> *invgf, int *iold, double *deltaold, int ipart_ref );
    //
    ////! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    //void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int bin ) override final;
    //
    ////! Project global current densities if Ionization in Species::dynamics,
    //void ionizationCurrents( Field *Jl, Field *Jr, Field *Jt, Particles &particles, int ipart, LocalFields Jion ) override final;
    //
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell,  int ipart_ref ) override final;
    //
    // Project susceptibility
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override final;
    
private:
};

#endif

