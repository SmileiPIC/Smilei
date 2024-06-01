#ifndef SMILEI_PROJECTOR_PROJECTOR1D2ORDERGPU_H
#define SMILEI_PROJECTOR_PROJECTOR1D2ORDERGPU_H

#include "Projector1D.h"


class Projector1D2OrderGPU : public Projector1D
{
public:
    Projector1D2OrderGPU( Params &parameters, Patch *a_patch );
    ~Projector1D2OrderGPU();

    /// For initialization and diags, doesn't use the standard scheme
    void basic( double      *rhoj,
                Particles   &particles,
                unsigned int ipart,
                unsigned int type,
                int bin_shift = 0 ) override;
    /// Projection wrapper
    void currentsAndDensityWrapper( ElectroMagn *EMfields,
                                    Particles   &particles,
                                    SmileiMPI   *smpi,
                                    int          istart,
                                    int          iend,
                                    int          ithread,
                                    bool         diag_flag,
                                    bool         is_spectral,
                                    int          ispec,
                                    int          icell     = 0,
                                    int          ipart_ref = 0 ) override;

    void susceptibility( ElectroMagn *EMfields,
                         Particles   &particles,
                         double       species_mass,
                         SmileiMPI   *smpi,
                         int          istart,
                         int          iend,
                         int          ithread,
                         int          icell     = 0,
                         int          ipart_ref = 0 ) override;
                         
    void ionizationCurrents( Field      *Jx,
                             Field      *Jy,
                             Field      *Jz,
                             Particles  &particles,
                             int         ipart,
                             LocalFields Jion ) override;


    //!Wrapper for task-based implementation of Smilei
    //! compiler complains otherwise even if it is completely useless
    void currentsAndDensityWrapperOnBuffers( double *b_Jx,
                                             double *b_Jy,
                                             double *b_Jz,
                                             double *b_rho,
                                             int bin_width,
                                             Particles &particles,
                                             SmileiMPI *smpi,
                                             int istart,
                                             int iend,
                                             int ithread,
                                             bool diag_flag,
                                             bool is_spectral,
                                             int ispec,
                                             int icell = 0,
                                             int ipart_ref = 0 ) override {};
/*#if defined( SMILEI_ACCELERATOR_MODE )

extern "C" void
currentDepositionKernel1DOnDevice( double *__restrict__ Jx,
                         double *__restrict__ Jy,
                         double *__restrict__ Jz,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         const double *__restrict__ particle_position_x,
                         const double *__restrict__ particle_momentum_y,
                         const double *__restrict__ particle_momentum_z,
                         const short *__restrict__ particle_charge,
                         const double *__restrict__ particle_weight,
                         const int *__restrict__ host_bin_index,
                         unsigned int x_dimension_bin_count,
                         const double *__restrict__ invgf_,
                         const int *__restrict__ iold_,
                         const double *__restrict__ deltaold_,
                         double inv_cell_volume,
                         double dx_inv,
                         double dx_ov_dt,
                         int    i_domain_begin,
                         int    not_spectral_ );

extern "C" void
currentAndDensityDepositionKernel1DOnDevice( double *__restrict__ Jx,
                                   double *__restrict__ Jy,
                                   double *__restrict__ Jz,
                                   double *__restrict__ rho,
                                   int Jx_size,
                                   int Jy_size,
                                   int Jz_size,
                                   int rho_size,
                                   const double *__restrict__ particle_position_x,
                                   const double *__restrict__ particle_momentum_y,
                                   const double *__restrict__ particle_momentum_z,
                                   const short *__restrict__ particle_charge,
                                   const double *__restrict__ particle_weight,
                                   const int *__restrict__ host_bin_index,
                                   unsigned int x_dimension_bin_count,
                                   const double *__restrict__ invgf_,
                                   const int *__restrict__ iold_,
                                   const double *__restrict__ deltaold_,
                                   double inv_cell_volume,
                                   double dx_inv,
                                   double dx_ov_dt,
                                   int    i_domain_begin,
                                   int    not_spectral_ );

#endif*/


protected:
    double dts2_;
    double dts4_;
    int    not_spectral_;
    unsigned int x_dimension_bin_count_;
};

#endif