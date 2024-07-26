//! HIP CUDA implementation

#ifndef Projector2D2OrderGPUKernelCUDAHIP_H
#define Projector2D2OrderGPUKernelCUDAHIP_H

#if defined( SMILEI_ACCELERATOR_GPU )

#if defined( __HIP__ )
    #include <hip/hip_runtime.h>
#elif defined( __NVCC__ )
    #include <cuda_runtime.h>
    #include <cuda.h>
#endif

#include "Params.h"
#include "gpu.h"



namespace cudahip2d {
//static
void currentDepositionKernel2D( double *__restrict__ host_Jx,
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
                               int    not_spectral_ );

//static 
void currentAndDensityDepositionKernel2D(
                                double *__restrict__ host_Jx,
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
                                int    not_spectral_ );

} // namespace cudahip2d

#endif
#endif

