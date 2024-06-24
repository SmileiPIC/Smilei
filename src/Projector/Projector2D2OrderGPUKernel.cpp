#if defined( SMILEI_ACCELERATOR_GPU )

#include "Projector2D2OrderGPUKernelCUDAHIP.h"
#include <cmath>
#include "Tools.h"
#include "gpu.h"

//! Project global current densities (EMfields->Jx_/Jy_/Jz_)
//!
extern "C" void
currentDepositionKernel2DOnDevice( double *__restrict__ host_Jx,
                         double *__restrict__ host_Jy,
                         double *__restrict__ host_Jz,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         const double *__restrict__ device_particle_position_x,
                         const double *__restrict__ device_particle_position_y,
                         const double *__restrict__ device_particle_momentum_z,
                         const short *__restrict__ device_particle_charge,
                         const double *__restrict__ device_particle_weight,
                         const int *__restrict__ host_bin_index,
                         unsigned int x_dimension_bin_count,
                         unsigned int y_dimension_bin_count,
                         const double *__restrict__ host_invgf_,
                         const int *__restrict__ host_iold_,
                         const double *__restrict__ host_deltaold_,
                         double inv_cell_volume,
                         double dx_inv,
                         double dy_inv,
                         double dx_ov_dt,
                         double dy_ov_dt,
                         int    i_domain_begin,
                         int    j_domain_begin,
                         int    nprimy,
                         int    not_spectral_ )
{
    //#if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    //naive:: // the naive, OMP version serves as a reference along with the CPU version
    //#else
    cudahip2d::
    //#endif
        currentDepositionKernel2D( host_Jx, host_Jy, host_Jz,
                                 Jx_size, Jy_size, Jz_size,
                                 device_particle_position_x, device_particle_position_y,
                                 device_particle_momentum_z,
                                 device_particle_charge,
                                 device_particle_weight,
                                 host_bin_index,
                                 x_dimension_bin_count,
                                 y_dimension_bin_count,
                                 host_invgf_,
                                 host_iold_, host_deltaold_,
                                 inv_cell_volume,
                                 dx_inv, dy_inv,
                                 dx_ov_dt, dy_ov_dt,
                                 i_domain_begin, j_domain_begin,
                                 nprimy,
                                 not_spectral_ );
}


//! Project global current and charge densities (EMfields->Jx_/Jy_/Jz_/rho_)
//!
extern "C" void
currentAndDensityDepositionKernel2DOnDevice( double *__restrict__ host_Jx,
                                   double *__restrict__ host_Jy,
                                   double *__restrict__ host_Jz,
                                   double *__restrict__ host_rho,
                                   int Jx_size,
                                   int Jy_size,
                                   int Jz_size,
                                   int rho_size,
                                   const double *__restrict__ device_particle_position_x,
                                   const double *__restrict__ device_particle_position_y,
                                   const double *__restrict__ device_particle_momentum_z,
                                   const short *__restrict__ device_particle_charge,
                                   const double *__restrict__ device_particle_weight,
                                   const int *__restrict__ host_bin_index,
                                   unsigned int x_dimension_bin_count,
                                   unsigned int y_dimension_bin_count,
                                   const double *__restrict__ host_invgf_,
                                   const int *__restrict__ host_iold_,
                                   const double *__restrict__ host_deltaold_,
                                   double inv_cell_volume,
                                   double dx_inv,
                                   double dy_inv,
                                   double dx_ov_dt,
                                   double dy_ov_dt,
                                   int    i_domain_begin,
                                   int    j_domain_begin,
                                   int    nprimy,
                                   int    not_spectral_ )
{
    //#if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    //naive:: // the naive, OMP version serves as a reference along with the CPU version
    //#else
    cudahip2d::
    //#endif
        currentAndDensityDepositionKernel2D( host_Jx, host_Jy, host_Jz, host_rho,
                                           Jx_size, Jy_size, Jz_size, rho_size,
                                           device_particle_position_x, device_particle_position_y,
                                           device_particle_momentum_z,
                                           device_particle_charge,
                                           device_particle_weight,
                                           host_bin_index,
                                           x_dimension_bin_count,
                                           y_dimension_bin_count,
                                           host_invgf_,
                                           host_iold_, host_deltaold_,
                                           inv_cell_volume,
                                           dx_inv, dy_inv,
                                           dx_ov_dt, dy_ov_dt,
                                           i_domain_begin, j_domain_begin,
                                           nprimy,
                                           not_spectral_ );
}
#endif

