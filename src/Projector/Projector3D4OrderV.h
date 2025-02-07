#ifndef PROJECTOR3D4ORDERV_H
#define PROJECTOR3D4ORDERV_H

#include "Projector3D.h"
#include "Pragma.h"

class Projector3D4OrderV : public Projector3D
{
public:

    //! Creator for Projector3D4OrderV
    Projector3D4OrderV( Params &, Patch *patch );

    ~Projector3D4OrderV();

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_)
    //! \param buffer_size number of particles in the buffers invgf, iold, deltaold
    inline void __attribute__((always_inline)) currents( double    * __restrict__ Jx,
                                                         double    * __restrict__ Jy,
                                                         double    * __restrict__ Jz,
                                                         Particles &particles,
                                                         unsigned int istart,
                                                         unsigned int iend,
                                                         double    * __restrict__ invgf,
                                                         int       * __restrict__ iold,
                                                         double    * __restrict__ deltaold,
                                                         unsigned int buffer_size,
                                                         int ipart_ref = 0, int bin_shift = 0 );

    //! Project global current densities (EMfields->Jx_/Jy_/Jz_/rho), diagFields timestep
    inline void __attribute__((always_inline)) currentsAndDensity( double *Jx,
                                                double  * __restrict__ Jy,
                                                double  * __restrict__ Jz,
                                                double  * __restrict__ rho,
                                                Particles &particles,
                                                unsigned int istart,
                                                unsigned int iend,
                                                double  * __restrict__ invgf,
                                                int     * __restrict__ iold,
                                                double  * __restrict__ deltaold,
                                                unsigned int buffer_size,
                                                int ipart_ref = 0, int bin_shift = 0 );

    //! Project global current charge (EMfields->rho_), frozen & diagFields timestep
    void basic( double *rhoj, Particles &particles, unsigned int ipart, unsigned int bin, int bin_shift = 0 ) override final;
    
    //! Project global current densities if Ionization in SpeciesV::dynamics,
    void ionizationCurrents( Field *Jx, Field *Jy, Field *Jz, Particles &particles, int ipart, LocalFields Jion ) override final;

    //!Wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, bool diag_flag, bool is_spectral, int ispec, int icell, int ipart_ref ) override final;
    
    void susceptibility( ElectroMagn *EMfields, Particles &particles, double species_mass, SmileiMPI *smpi, int istart, int iend,  int ithread, int icell, int ipart_ref ) override;

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
