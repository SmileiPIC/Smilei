#ifndef PROJECTOR2D2ORDERV_H
#define PROJECTOR2D2ORDERV_H

#include "Projector2D.h"

class Projector2D2OrderV : public Projector2D
{
public:

    // Creator
    Projector2D2OrderV( Params &, Patch *patch );

    // Destructor
    ~Projector2D2OrderV();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! \param buffer_size number of particles in the buffers invgf, iold, deltaold
    void currents( double * __restrict__ Jx,
                   double * __restrict__ Jy,
                   double * __restrict__ Jz,
                   Particles &particles,
                   unsigned int istart, unsigned int iend,
                   double * __restrict__ invgf,
                   int    * __restrict__ iold,
                   double * __restrict__ deltaold,
                   unsigned int buffer_size, int ipart_ref = 0 );

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    //! \param buffer_size number of particles in the buffers invgf, iold, deltaold
    inline void __attribute__((always_inline)) currentsAndDensity(  double * __restrict__ Jx,
                                                                    double * __restrict__ Jy,
                                                                    double * __restrict__ Jz,
                                                                    double * __restrict__ rho,
                                                                    Particles &particles, unsigned int istart,
                                                                    unsigned int iend,
                                                                    double * __restrict__ invgf,
                                                                    int    * __restrict__ iold,
                                                                    double * __restrict__ deltaold,
                                                                    unsigned int buffer_size,
                                                                    int ipart_ref );

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int bin ) override final;

    //! Project global current densities if Ionization in Species::dynamics,
    void ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion ) override final;

    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell,  int ipart_ref ) override final;

    // Project susceptibility
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override final;

private:
};

#endif
