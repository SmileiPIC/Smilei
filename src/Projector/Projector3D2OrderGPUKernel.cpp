
// TODO(Etienne M): The makefile does not recognise this file and doesn't compute
// it's dependencies. If you make a modification in one of the header this file
// includes, you must `touch` this file. IF you dont do that you'll have ABI/ODR
// issues (!).

#if defined( SMILEI_ACCELERATOR_MODE )

    //! Simple switch to jump between the reference (omp) implementation and the
    //! hip one.
    //! NOTE: If you wanna use the OMP version, you must rename this file to
    //! .cpp instead of .cu for the HIP. The preprocessor and the Smilei
    //! makefile will take care of the rest.
    //!
    #if defined( __HIP__ ) // || defined (__CUDACC__)
    // HIP compiler support enabled (for .cu files)
    #else
        #define PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION 1
    #endif

    // #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    //     #include <cmath>

    //     #include "Tools.h"
    //     #include "gpu.h"
    // // #elif defined( __CUDACC__ ) 
    // //     #include "Params.h"
    // //     #include "gpu.h"
    // #elif defined( __HIP__ )
    //     #include <hip/hip_runtime.h>

    //     #include "Params.h"
    //     #include "gpu.h"
    // #endif

    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )

    // #include "Projector3D2OrderGPUKernelAcc.h"
    #include "Projector3D2OrderGPUKernelNaive.h"

    #else

    #include "Projector3D2OrderGPUKernelHIP.h"

    #endif

//! Project global current densities (EMfields->Jx_/Jy_/Jz_)
//!
extern "C" void
currentDepositionKernel3D( double *__restrict__ host_Jx,
                           double *__restrict__ host_Jy,
                           double *__restrict__ host_Jz,
                           int Jx_size,
                           int Jy_size,
                           int Jz_size,
                           const double *__restrict__ device_particle_position_x,
                           const double *__restrict__ device_particle_position_y,
                           const double *__restrict__ device_particle_position_z,
                           const short  *__restrict__ device_particle_charge,
                           const double *__restrict__ device_particle_weight,
                           const int    *__restrict__ host_bin_index,
                           unsigned int x_dimension_bin_count,
                           unsigned int y_dimension_bin_count,
                           unsigned int z_dimension_bin_count,
                           const double *__restrict__ host_invgf_,
                           int    *__restrict__ host_iold,
                           const double *__restrict__ host_deltaold_,
                           const unsigned int number_of_particles,
                           double inv_cell_volume,
                           double dx_inv,
                           double dy_inv,
                           double dz_inv,
                           double dx_ov_dt,
                           double dy_ov_dt,
                           double dz_ov_dt,
                           int    i_domain_begin,
                           int    j_domain_begin,
                           int    k_domain_begin,
                           int    nprimy,
                           int    nprimz,
                           int    not_spectral )
{
    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    acc:: // the naive, OMP version serves as a reference along with the CPU version
    #else
    hip::
    #endif
        currentDepositionKernel3D( host_Jx, host_Jy, host_Jz,
                                   Jx_size, Jy_size, Jz_size,
                                   device_particle_position_x, device_particle_position_y,
                                   device_particle_position_z,
                                   device_particle_charge,
                                   device_particle_weight,
                                   host_bin_index,
                                   x_dimension_bin_count,
                                   y_dimension_bin_count,
                                   z_dimension_bin_count,
                                   host_invgf_,
                                   host_iold,
                                   host_deltaold_,
                                   number_of_particles,
                                   inv_cell_volume,
                                   dx_inv, dy_inv, dz_inv,
                                   dx_ov_dt, dy_ov_dt, dz_ov_dt,
                                   i_domain_begin, j_domain_begin, k_domain_begin,
                                   nprimy, nprimz,
                                   not_spectral );
}

//! Project global current and charge densities (EMfields->Jx_/Jy_/Jz_/rho_)
//!
extern "C" void
currentAndDensityDepositionKernel3D( double *__restrict__ host_Jx,
                                     double *__restrict__ host_Jy,
                                     double *__restrict__ host_Jz,
                                     double *__restrict__ host_rho,
                                     int Jx_size,
                                     int Jy_size,
                                     int Jz_size,
                                     int rho_size,
                                     const double *__restrict__ device_particle_position_x,
                                     const double *__restrict__ device_particle_position_y,
                                     const double *__restrict__ device_particle_position_z,
                                     const short *__restrict__ device_particle_charge,
                                     const double *__restrict__ device_particle_weight,
                                     const int *__restrict__ host_bin_index,
                                     unsigned int x_dimension_bin_count,
                                     unsigned int y_dimension_bin_count,
                                     unsigned int z_dimension_bin_count,
                                     const double *__restrict__ host_invgf_,
                                     int *__restrict__ host_iold,
                                     const double *__restrict__ host_deltaold_,
                                     const unsigned int number_of_particles,
                                     double inv_cell_volume,
                                     double dx_inv,
                                     double dy_inv,
                                     double dz_inv,
                                     double dx_ov_dt,
                                     double dy_ov_dt,
                                     double dz_ov_dt,
                                     int    i_domain_begin,
                                     int    j_domain_begin,
                                     int    k_domain_begin,
                                     int    nprimy,
                                     int    nprimz,
                                     int    not_spectral )
{
    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    acc:: // the naive, OMP version serves as a reference along with the CPU version
    #else
    hip::
    #endif
        currentAndDensityDepositionKernel3D( host_Jx, host_Jy, host_Jz, host_rho,
                                             Jx_size, Jy_size, Jz_size, rho_size,
                                             device_particle_position_x, 
                                             device_particle_position_y,
                                             device_particle_position_z,
                                             device_particle_charge,
                                             device_particle_weight,
                                             host_bin_index,
                                             x_dimension_bin_count,
                                             y_dimension_bin_count,
                                             z_dimension_bin_count,
                                             host_invgf_,
                                             host_iold,
                                             host_deltaold_,
                                             number_of_particles,
                                             inv_cell_volume,
                                             dx_inv, dy_inv, dz_inv,
                                             dx_ov_dt, dy_ov_dt, dz_ov_dt,
                                             i_domain_begin, j_domain_begin, k_domain_begin,
                                             nprimy, nprimz,
                                             not_spectral );
}

#endif
