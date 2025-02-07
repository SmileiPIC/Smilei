#ifndef PROJECTOR3D2ORDER_H
#define PROJECTOR3D2ORDER_H

#include "Projector3D.h"


class Projector3D2Order : public Projector3D
{
public:
    Projector3D2Order( Params &, Patch *patch );
    ~Projector3D2Order();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    inline void __attribute__((always_inline)) currents( double *Jx, double *Jy, double *Jz, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, int bin_shift = 0 );
    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void __attribute__((always_inline)) currentsAndDensity( double *Jx, double *Jy, double *Jz, double *rho, Particles &particles, unsigned int ipart, double invgf, int *iold, double *deltaold, int bin_shift = 0 );
    
    //! Project global current charge (EMfields->rho_ , J), for initialization and diags
    void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int type, int bin_shift = 0 ) override final;
    
    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion ) override final;
    
    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell = 0, int ipart_ref = 0 ) override final;
    
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell = 0, int ipart_ref = 0 ) override final;

private:
    int pxr;
    double dt, dts2, dts4;
};

#endif
