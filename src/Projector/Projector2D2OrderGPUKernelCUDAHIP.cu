
// TODO(Etienne M): The makefile does not recognise this file and doesn't compute
// it's dependencies. If you make a modification in one of the header this file
// includes, you must `touch` this file. IF you dont do that you'll have ABI/ODR
// issues (!).

//#if defined( SMILEI_ACCELERATOR_GPU_OMP )


#if defined( __HIP__ ) 
    #include <hip/hip_runtime.h>
#elif defined( __NVCC__ )
    #include <cuda_runtime.h>
    #include <cuda.h>
#endif

#include "Params.h"
#include "gpu.h"
#include <iostream>

#if defined( __HIP__ )
  // HIP compiler support enabled (for .cu files)
#else
    #define PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION 1
#endif

#if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
    #include <cmath>

    #include "Tools.h"
#else
    #include <hip/hip_runtime.h>

    #include "Params.h"
    #include "gpu.h"
#endif

//    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )

//namespace naive {
//
//    void //static inline void
//    currentDepositionKernel2D( double     *__restrict__ Jx,
//                             double       *__restrict__ Jy,
//                             double       *__restrict__ Jz,
//                             int Jx_size,
//                             int Jy_size,
//                             int Jz_size,
//                             const double *__restrict__ device_particle_position_x,
//                             const double *__restrict__ device_particle_position_y,
//                             const double *__restrict__ device_particle_momentum_z,
//                             const short  *__restrict__ device_particle_charge,
//                             const double *__restrict__ device_particle_weight,
//                             const int    *__restrict__ host_bin_index,
//                             unsigned int x_dimension_bin_count,
//                             unsigned int y_dimension_bin_count,
//                             const double *__restrict__ invgf_,
//                             const int *__restrict__ iold_,
//                             const double *__restrict__ deltaold_,
//                             double inv_cell_volume,
//                             double dx_inv,
//                             double dy_inv,
//                             double dx_ov_dt,
//                             double dy_ov_dt,
//                             int    i_domain_begin,
//                             int    j_domain_begin,
//                             int    nprimy,
//                             int    not_spectral_ )
//    {
//        // The OMP implementation is NOT bin aware. As per the precondition on
//        // host_bin_index, index zero always contains the number of particles.
//        // See nvidiaParticles::prepareBinIndex / setHostBinIndex.
//        const unsigned int bin_count      = 1;
//        const int          particle_count = host_bin_index[bin_count - 1];
//
//        #if defined( SMILEI_ACCELERATOR_GPU_OMP )
//            #pragma omp target is_device_ptr /* map */ ( /* to: */                                            \
//                                                         device_particle_position_x /* [0:particle_count] */, \
//                                                         device_particle_position_y /* [0:particle_count] */, \
//                                                         device_particle_momentum_z /* [0:particle_count] */, \
//                                                         device_particle_charge /* [0:particle_count] */,     \
//                                                         device_particle_weight /* [0:particle_count] */ )
//            #pragma omp teams thread_limit( 64 ) distribute parallel for
//        #elif defined( SMILEI_ACCELERATOR_GPU_OACC )
//            #pragma acc parallel                      \
//            deviceptr( device_particle_position_x,    \
//                       device_particle_position_y,    \
//                       device_particle_momentum_z,    \
//                       device_particle_charge,        \
//                       device_particle_weight )       \
//                present( iold [0:3 * particle_count], \
//                         deltaold [0:3 * particle_count] )
//            #pragma acc loop gang worker vector
//        #endif
//        for( int particle_index = 0; particle_index < particle_count; ++particle_index ) {
//            const double invgf                        = invgf_[particle_index];
//            const int *const __restrict__ iold        = &iold_[particle_index];
//            const double *const __restrict__ deltaold = &deltaold_[particle_index];
//
//            double Sx0[5];
//            double Sx1[5];
//            double Sy0[5];
//            double Sy1[5];
//
//            // Variable declaration & initialization
//            // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf
//
//            // Locate the particle on the primal grid at former time-step & calculate coeff. S0
//            {
//                const double delta  = deltaold[0 * particle_count];
//                const double delta2 = delta * delta;
//                Sx0[0]              = 0.0;
//                Sx0[1]              = 0.5 * ( delta2 - delta + 0.25 );
//                Sx0[2]              = 0.75 - delta2;
//                Sx0[3]              = 0.5 * ( delta2 + delta + 0.25 );
//                Sx0[4]              = 0.0;
//            }
//            {
//                const double delta  = deltaold[1 * particle_count];
//                const double delta2 = delta * delta;
//                Sy0[0]              = 0.0;
//                Sy0[1]              = 0.5 * ( delta2 - delta + 0.25 );
//                Sy0[2]              = 0.75 - delta2;
//                Sy0[3]              = 0.5 * ( delta2 + delta + 0.25 );
//                Sy0[4]              = 0.0;
//            }
//
//            // Locate the particle on the primal grid at current time-step & calculate coeff. S1
//            {
//                const double xpn      = device_particle_position_x[particle_index] * dx_inv;
//                const int    ip       = std::round( xpn );
//                const int    ipo      = iold[0 * particle_count];
//                const int    ip_m_ipo = ip - ipo - i_domain_begin;
//                const double delta    = xpn - static_cast<double>( ip );
//                const double delta2   = delta * delta;
//
//                Sx1[0] = 0.0;
//                Sx1[1] = 0.0;
//                // Sx1[2] = 0.0; // Always set below
//                Sx1[3] = 0.0;
//                Sx1[4] = 0.0;
//
//                Sx1[ip_m_ipo + 1] = 0.5 * ( delta2 - delta + 0.25 );
//                Sx1[ip_m_ipo + 2] = 0.75 - delta2;
//                Sx1[ip_m_ipo + 3] = 0.5 * ( delta2 + delta + 0.25 );
//            }
//            {
//                const double ypn      = device_particle_position_y[particle_index] * dy_inv;
//                const int    jp       = std::round( ypn );
//                const int    jpo      = iold[1 * particle_count];
//                const int    jp_m_jpo = jp - jpo - j_domain_begin;
//                const double delta    = ypn - static_cast<double>( jp );
//                const double delta2   = delta * delta;
//
//                Sy1[0] = 0.0;
//                Sy1[1] = 0.0;
//                // Sy1[2] = 0.0; // Always set below
//                Sy1[3] = 0.0;
//                Sy1[4] = 0.0;
//
//                Sy1[jp_m_jpo + 1] = 0.5 * ( delta2 - delta + 0.25 );
//                Sy1[jp_m_jpo + 2] = 0.75 - delta2;
//                Sy1[jp_m_jpo + 3] = 0.5 * ( delta2 + delta + 0.25 );
//            }
//
//            // (x,y,z) components of the current density for the macro-particle
//            const double charge_weight = inv_cell_volume * static_cast<double>( device_particle_charge[particle_index] ) * device_particle_weight[particle_index];
//            const double crx_p         = charge_weight * dx_ov_dt;
//            const double cry_p         = charge_weight * dy_ov_dt;
//            const double crz_p         = charge_weight * ( 1.0 / 3.0 ) * device_particle_momentum_z[particle_index] * invgf;
//
//            // This is the particle position as grid index
//            // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
//            const int ipo = iold[0 * particle_count] - 2;
//            const int jpo = iold[1 * particle_count] - 2;
//
//            for( unsigned int i = 0; i < 1; ++i ) {
//                const int iloc = ( i + ipo ) * nprimy + jpo;
//                /* Jx[iloc] += tmpJx[0]; */
//
//                SMILEI_ACCELERATOR_ATOMIC
//                Jz[iloc] += crz_p * ( Sy1[0] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
//                double tmp = 0.0;
//                for( unsigned int j = 1; j < 5; j++ ) {
//                    tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jy[iloc + j + not_spectral_ * ( /* i + */ ipo )] += tmp;
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] /* + Sx0[i] */ ) +
//                                              Sy1[j] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
//                }
//            }
//
//            double tmpJx[5]{};
//
//            for( unsigned int i = 1; i < 5; ++i ) {
//                const int iloc = ( i + ipo ) * nprimy + jpo;
//                tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( 0.5 * ( Sy1[0] - Sy0[0] ) );
//                SMILEI_ACCELERATOR_ATOMIC
//                Jx[iloc] += tmpJx[0];
//                SMILEI_ACCELERATOR_ATOMIC
//                Jz[iloc] += crz_p * ( Sy1[0] * ( 0.5 * Sx0[i] + Sx1[i] ) );
//                double tmp = 0.0;
//                for( unsigned int j = 1; j < 5; ++j ) {
//                    tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + 0.5 * ( Sy1[j] - Sy0[j] ) );
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jx[iloc + j] += tmpJx[j];
//                    tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jy[iloc + j + not_spectral_ * ( i + ipo )] += tmp;
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] + Sx0[i] ) +
//                                              Sy1[j] * ( 0.5 * Sx0[i] + Sx1[i] ) );
//                }
//            }
//        }
//    } // end currentDepositionKernel
//
//    //static inline
//    void
//    currentAndDensityDepositionKernel( double *__restrict__ Jx,
//                                       double *__restrict__ Jy,
//                                       double *__restrict__ Jz,
//                                       double *__restrict__ rho,
//                                       int Jx_size,
//                                       int Jy_size,
//                                       int Jz_size,
//                                       int rho_size,
//                                       const double *__restrict__ device_particle_position_x,
//                                       const double *__restrict__ device_particle_position_y,
//                                       const double *__restrict__ device_particle_momentum_z,
//                                       const short *__restrict__ device_particle_charge,
//                                       const double *__restrict__ device_particle_weight,
//                                       const int *__restrict__ host_bin_index,
//                                       unsigned int,
//                                       unsigned int,
//                                       const double *__restrict__ invgf_,
//                                       const int *__restrict__ iold_,
//                                       const double *__restrict__ deltaold_,
//                                       double inv_cell_volume,
//                                       double dx_inv,
//                                       double dy_inv,
//                                       double dx_ov_dt,
//                                       double dy_ov_dt,
//                                       int    i_domain_begin,
//                                       int    j_domain_begin,
//                                       int    nprimy,
//                                       int    not_spectral_ )
//    {
//        // The OMP implementation is NOT bin aware. As per the precondition on
//        // host_bin_index, index zero always contains the number of particles.
//        // See nvidiaParticles::prepareBinIndex / setHostBinIndex.
//        const unsigned int bin_count      = 1;
//        const int          particle_count = host_bin_index[bin_count - 1];
//
//        #if defined( SMILEI_ACCELERATOR_GPU_OMP )
//            #pragma omp target is_device_ptr /* map */ ( /* to: */                                            \
//                                                         device_particle_position_x /* [0:particle_count] */, \
//                                                         device_particle_position_y /* [0:particle_count] */, \
//                                                         device_particle_momentum_z /* [0:particle_count] */, \
//                                                         device_particle_charge /* [0:particle_count] */,     \
//                                                         device_particle_weight /* [0:particle_count] */ )
//            #pragma omp teams thread_limit( 64 ) distribute parallel for
//        #elif defined( SMILEI_ACCELERATOR_GPU_OACC )
//            #pragma acc parallel                      \
//            deviceptr( device_particle_position_x,    \
//                       device_particle_position_y,    \
//                       device_particle_momentum_z,    \
//                       device_particle_charge,        \
//                       device_particle_weight )       \
//                present( iold [0:3 * particle_count], \
//                         deltaold [0:3 * particle_count] )
//            #pragma acc loop gang worker vector
//        #endif
//        for( int particle_index = 0; particle_index < particle_count; ++particle_index ) {
//            const double invgf                        = invgf_[particle_index];
//            const int *const __restrict__ iold        = &iold_[particle_index];
//            const double *const __restrict__ deltaold = &deltaold_[particle_index];
//
//            double Sx0[5];
//            double Sx1[5];
//            double Sy0[5];
//            double Sy1[5];
//
//            // Variable declaration & initialization
//            // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf
//
//            // Locate the particle on the primal grid at former time-step & calculate coeff. S0
//            {
//                const double delta  = deltaold[0 * particle_count];
//                const double delta2 = delta * delta;
//                Sx0[0]              = 0.0;
//                Sx0[1]              = 0.5 * ( delta2 - delta + 0.25 );
//                Sx0[2]              = 0.75 - delta2;
//                Sx0[3]              = 0.5 * ( delta2 + delta + 0.25 );
//                Sx0[4]              = 0.0;
//            }
//            {
//                const double delta  = deltaold[1 * particle_count];
//                const double delta2 = delta * delta;
//                Sy0[0]              = 0.0;
//                Sy0[1]              = 0.5 * ( delta2 - delta + 0.25 );
//                Sy0[2]              = 0.75 - delta2;
//                Sy0[3]              = 0.5 * ( delta2 + delta + 0.25 );
//                Sy0[4]              = 0.0;
//            }
//
//            // Locate the particle on the primal grid at current time-step & calculate coeff. S1
//            {
//                const double xpn      = device_particle_position_x[particle_index] * dx_inv;
//                const int    ip       = std::round( xpn );
//                const int    ipo      = iold[0 * particle_count];
//                const int    ip_m_ipo = ip - ipo - i_domain_begin;
//                const double delta    = xpn - static_cast<double>( ip );
//                const double delta2   = delta * delta;
//
//                Sx1[0] = 0.0;
//                Sx1[1] = 0.0;
//                // Sx1[2] = 0.0; // Always set below
//                Sx1[3] = 0.0;
//                Sx1[4] = 0.0;
//
//                Sx1[ip_m_ipo + 1] = 0.5 * ( delta2 - delta + 0.25 );
//                Sx1[ip_m_ipo + 2] = 0.75 - delta2;
//                Sx1[ip_m_ipo + 3] = 0.5 * ( delta2 + delta + 0.25 );
//            }
//            {
//                const double ypn      = device_particle_position_y[particle_index] * dy_inv;
//                const int    jp       = std::round( ypn );
//                const int    jpo      = iold[1 * particle_count];
//                const int    jp_m_jpo = jp - jpo - j_domain_begin;
//                const double delta    = ypn - static_cast<double>( jp );
//                const double delta2   = delta * delta;
//
//                Sy1[0] = 0.0;
//                Sy1[1] = 0.0;
//                // Sy1[2] = 0.0; // Always set below
//                Sy1[3] = 0.0;
//                Sy1[4] = 0.0;
//
//                Sy1[jp_m_jpo + 1] = 0.5 * ( delta2 - delta + 0.25 );
//                Sy1[jp_m_jpo + 2] = 0.75 - delta2;
//                Sy1[jp_m_jpo + 3] = 0.5 * ( delta2 + delta + 0.25 );
//            }
//
//            // (x,y,z) components of the current density for the macro-particle
//            const double charge_weight = inv_cell_volume * static_cast<double>( device_particle_charge[particle_index] ) * device_particle_weight[particle_index];
//            const double crx_p         = charge_weight * dx_ov_dt;
//            const double cry_p         = charge_weight * dy_ov_dt;
//            const double crz_p         = charge_weight * ( 1.0 / 3.0 ) * device_particle_momentum_z[particle_index] * invgf;
//
//            // This is the particle position as grid index
//            // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
//            const int ipo = iold[0 * particle_count] - 2;
//            const int jpo = iold[1 * particle_count] - 2;
//
//            // case i =0
//            for( unsigned int i = 0; i < 1; ++i ) {
//                const int iloc = ( i + ipo ) * nprimy + jpo;
//                /* Jx[iloc] += tmpJx[0]; */
//
//                SMILEI_ACCELERATOR_ATOMIC
//                Jz[iloc] += crz_p * ( Sy1[0] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
//
//                SMILEI_ACCELERATOR_ATOMIC
//                rho[iloc] += charge_weight * Sx1[0] * Sy1[0];
//                double tmp = 0.0;
//                for( unsigned int j = 1; j < 5; j++ ) {
//                    tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jy[iloc + j + not_spectral_ * ( /* i + */ ipo )] += tmp;
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] /* + Sx0[i] */ ) +
//                                              Sy1[j] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
//                    SMILEI_ACCELERATOR_ATOMIC
//                    rho[iloc + j] += charge_weight * Sx1[0] * Sy1[j];
//                }
//            }
//
//            double tmpJx[5]{};
//
//            // case i> 0
//            for( unsigned int i = 1; i < 5; ++i ) {
//                const int iloc = ( i + ipo ) * nprimy + jpo;
//                tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( 0.5 * ( Sy1[0] - Sy0[0] ) );
//
//                SMILEI_ACCELERATOR_ATOMIC
//                Jx[iloc] += tmpJx[0];
//
//                SMILEI_ACCELERATOR_ATOMIC
//                Jz[iloc] += crz_p * ( Sy1[0] * ( 0.5 * Sx0[i] + Sx1[i] ) );
//
//                SMILEI_ACCELERATOR_ATOMIC
//                rho[iloc] += charge_weight * Sx1[i] * Sy1[0];
//
//                double tmp = 0.0;
//                for( unsigned int j = 1; j < 5; ++j ) {
//                    tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + 0.5 * ( Sy1[j] - Sy0[j] ) );
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jx[iloc + j] += tmpJx[j];
//                    tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jy[iloc + j + not_spectral_ * ( i + ipo )] += tmp;
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] + Sx0[i] ) +
//                                              Sy1[j] * ( 0.5 * Sx0[i] + Sx1[i] ) );
//
//                    SMILEI_ACCELERATOR_ATOMIC
//                    rho[iloc + j] += charge_weight * Sx1[i] * Sy1[j];
//                }
//            }
//        }
//    } // end currentDepositionKernel
//
//
//} // namespace naive
//
//    #else

namespace cudahip2d {
    namespace detail {
#if defined( __HIP__ )
        static inline void
        checkErrors( ::hipError_t an_error_code,
                     const char  *file_name,
                     int          line )
        {
            if( an_error_code != ::hipError_t::hipSuccess ) {
                std::cout << "HIP error at " << file_name << ":" << line
                          << " -> " << ::hipGetErrorString( an_error_code ) << std::endl;
                std::exit( EXIT_FAILURE );
            }
        }
// For NVIDIA compiler 
#elif defined(  __NVCC__ )
        static inline void
        checkErrors( ::cudaError_t an_error_code,
                     const char  *file_name,
                     int          line )
        {
            if( an_error_code != ::cudaError_t::cudaSuccess ) {
                std::cout << "CUDA error at " << file_name << ":" << line << " -> " << ::cudaGetErrorString( an_error_code ) << std::endl;
                std::exit( EXIT_FAILURE );
            }
        }
#endif

   } // namespace detail

    #define checkHIPErrors( an_expression )                           \
        do {                                                          \
            detail::checkErrors( an_expression, __FILE__, __LINE__ ); \
        } while( 0 )  

    namespace kernel {
        namespace atomic {
            namespace LDS {
                __device__ void
                AddNoReturn( float *a_pointer, float a_value )
                {
        #if defined( __gfx90a__ )
                    ::unsafeAtomicAdd( a_pointer, a_value );

                    // uint32_t *as_uint32{ reinterpret_cast<uint32_t *>( a_pointer ) };
                    // uint32_t  last_seen_value{ __atomic_load_n( as_uint32, __ATOMIC_RELAXED ) };
                    // uint32_t  assumed;
                    // do {
                    //     assumed         = last_seen_value;
                    //     last_seen_value = ::atomicCAS( as_uint32,
                    //                                    last_seen_value,
                    //                                    __float_as_uint( a_value + __float_as_uint( last_seen_value ) ) );
                    // } while( assumed != last_seen_value );
        #else
                    ::atomicAdd( a_pointer, a_value );
        #endif
                }

                __device__ void
                AddNoReturn( double *a_pointer, double a_value )
                {
        #if defined( __gfx90a__ )
                    ::unsafeAtomicAdd( a_pointer, a_value );
        #else
                    ::atomicAdd( a_pointer, a_value );
        #endif
                }
            } // namespace LDS

            namespace GDS {
                __device__ void
                AddNoReturn( double *a_pointer, double a_value )
                {
        #if defined( __gfx90a__ )
                    ::unsafeAtomicAdd( a_pointer, a_value );
        #else
                    ::atomicAdd( a_pointer, a_value );
        #endif
                }
            } // namespace GDS
        }     // namespace atomic


        template <typename ComputeFloat>
        __device__ void inline __attribute__((always_inline)) init_S0(const ComputeFloat delta, ComputeFloat *__restrict__ S0)
        {
            const ComputeFloat delta2 = delta * delta;
            S0[0] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
            S0[1] = static_cast<ComputeFloat>( 0.75 ) - delta2;
            S0[2] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
            S0[3] = static_cast<ComputeFloat>( 0.0 ) ;
        }

        template <typename ComputeFloat>
        __device__ void inline __attribute__((always_inline)) init_S1(const ComputeFloat xpn, const int ipo,  const int i_domain_begin,
                                                                      ComputeFloat *__restrict__ S1)
        {
            // const int    ip        = static_cast<int>( xpn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
            const int          ip       = std::round( xpn );
            const int          ip_m_ipo = ip - ipo - i_domain_begin;
            const ComputeFloat delta    = xpn - static_cast<ComputeFloat>( ip );
            const ComputeFloat delta2   = delta * delta;

            S1[0] = static_cast<ComputeFloat>( 0.0 );
            S1[1] = static_cast<ComputeFloat>( 0.0 ); // S1[2] = 0.0; // Always set below
            S1[3] = static_cast<ComputeFloat>( 0.0 );
            S1[4] = static_cast<ComputeFloat>( 0.0 );

            S1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
            S1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
            S1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
        }


        template <typename ComputeFloat,
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize>
        __global__ void
        // __launch_bounds__(kWorkgroupSize, 1)
        DepositCurrentDensity_2D_Order2( double *__restrict__ device_Jx,
                                         double *__restrict__ device_Jy,
                                         double *__restrict__ device_Jz,
                                         int Jx_size,
                                         int Jy_size,
                                         int Jz_size,
                                         const double *__restrict__ device_particle_position_x,
                                         const double *__restrict__ device_particle_position_y,
                                         const double *__restrict__ device_particle_momentum_z,
                                         const short *__restrict__ device_particle_charge,
                                         const double *__restrict__ device_particle_weight,
                                         const int *__restrict__ device_bin_index,
                                         const double *__restrict__ device_invgf_,
                                         const int *__restrict__ device_iold_,
                                         const double *__restrict__ device_deltaold_,
                                         ComputeFloat inv_cell_volume,
                                         ComputeFloat dx_inv,
                                         ComputeFloat dy_inv,
                                         ComputeFloat dx_ov_dt,
                                         ComputeFloat dy_ov_dt,
                                         int          i_domain_begin,
                                         int          j_domain_begin,
                                         int          nprimy,
                                         int          not_spectral_ )
        {
            // TODO(Etienne M): refactor this function. Break it into smaller
            // pieces (lds init/store, coeff computation, deposition etc..)
            // TODO(Etienne M): __ldg could be used to slightly improve GDS load
            // speed. This would only have an effect on Nvidia cards as this
            // operation is a no op on AMD.
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x * gridDim.y;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y + y_cluster_coordinate; // The indexing order is: x * ywidth * zwidth + y * zwidth + z
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 2 /* 2D */ );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 2 /* 2D */ );
            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ );
            ComputeFloat one_half = 1. / 2.;

            // NOTE: We gain from the particles not being sorted inside a
            // cluster because it reduces the bank conflicts one gets when
            // multiple threads access the same part of the shared memory. Such
            // "conflicted" accesses are serialized !
            // NOTE: We use a bit to much LDS. For Jx, the first row could be
            // discarded, for Jy we could remove the first column.

            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 2 /* 2D */, 2 /* 2nd order interpolation */ );

            // NOTE: I tried having only one cache and reusing it. Doing that
            // requires you to iterate multiple time over the particle which is
            // possible but cost more bandwidth. The speedup was ~x0.92.
            __shared__ ReductionFloat Jx_scratch_space[kFieldScratchSpaceSize];
            __shared__ ReductionFloat Jy_scratch_space[kFieldScratchSpaceSize];
            __shared__ ReductionFloat Jz_scratch_space[kFieldScratchSpaceSize];

            // Init the shared memory

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {
                Jx_scratch_space[field_index] = static_cast<ReductionFloat>( 0.0 );
                Jy_scratch_space[field_index] = static_cast<ReductionFloat>( 0.0 );
                Jz_scratch_space[field_index] = static_cast<ReductionFloat>( 0.0 );
            }

            __syncthreads();

            const unsigned int particle_count = device_bin_index[bin_count - 1];

            // This workgroup has to process distance(last_particle,
            // first_particle) particles
            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 :
                                                                                     device_bin_index[workgroup_dedicated_bin_index - 1];
            const unsigned int last_particle  = device_bin_index[workgroup_dedicated_bin_index];

            for( unsigned int particle_index = first_particle + thread_index_offset;
                 particle_index < last_particle;
                 particle_index += loop_stride ) {
                const ComputeFloat invgf                  = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__ iold        = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];

                // Variable declaration & initialization
                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // Locate the particle on the primal grid at former time-step & calculate coeff. S0
                {
                    const ComputeFloat delta  = deltaold[0 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sx0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[4] = static_cast<ComputeFloat>( 0.0 );
                }
                {
                    const ComputeFloat delta  = deltaold[1 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sy0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[4] = static_cast<ComputeFloat>( 0.0 );
                }//*/
                //init_S0(deltaold[0 * particle_count], Sx0);
                //init_S0(deltaold[1 * particle_count], Sy0);


                // Locate the particle on the primal grid at current time-step & calculate coeff. S1
                {
                    // const int    ip             = static_cast<int>( xpn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat xpn      = static_cast<ComputeFloat>( device_particle_position_x[particle_index] ) * dx_inv;
                    const int          ip       = std::round( xpn );
                    const int          ipo      = iold[0 * particle_count];
                    const int          ip_m_ipo = ip - ipo - i_domain_begin;
                    const ComputeFloat delta    = xpn - static_cast<ComputeFloat>( ip );
                    const ComputeFloat delta2   = delta * delta;

                    Sx1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sx1[2] = 0.0; // Always set below
                    Sx1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sx1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }
                {
                    // const int    jp             = static_cast<int>( ypn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat ypn      = static_cast<ComputeFloat>( device_particle_position_y[particle_index] ) * dy_inv;
                    const int          jp       = std::round( ypn );
                    const int          jpo      = iold[1 * particle_count];
                    const int          jp_m_jpo = jp - jpo - j_domain_begin;
                    const ComputeFloat delta    = ypn - static_cast<ComputeFloat>( jp );
                    const ComputeFloat delta2   = delta * delta;

                    Sy1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sy1[2] = 0.0; // Always set below
                    Sy1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sy1[jp_m_jpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy1[jp_m_jpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy1[jp_m_jpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }//*/
                /*init_S1( static_cast<ComputeFloat>( device_particle_position_x[particle_index] ) * dx_inv,
                    iold[0 * particle_count], i_domain_begin, Sx1);
                init_S1( static_cast<ComputeFloat>( device_particle_position_y[particle_index] ) * dy_inv,
                    iold[1 * particle_count], j_domain_begin, Sy1);//*/


                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( 1.0 / 3.0 ) * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;

                // Jx

                ComputeFloat tmpJx[5]{};

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 , 2 ) + jpo;
                    tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( static_cast<ComputeFloat>( 0.5 ) * ( Sy1[0] - Sy0[0] ) );
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmpJx[0] ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + static_cast<ComputeFloat>( 0.5 ) * ( Sy1[j] - Sy0[j] ) );
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc + j], static_cast<ReductionFloat>( tmpJx[j] ) );
                    }
                }//*/

                /*{
                    ComputeFloat tmp = crx_p * static_cast<ComputeFloat>( 0.5 ) * Sy1[0] ;
                    ComputeFloat tmp_reduction{};
                    const int iloc = ipo *  GPUClusterWithGCWidth + jpo ;
                    tmp_reduction -= Sx1[0] * tmp;
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int i = 2; i < 5; ++i ) {
                        tmp_reduction -= ( Sx1[i - 1] - Sx0[i - 2] ) * tmp;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc + i * GPUClusterWithGCWidth], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }

                for( unsigned int j = 1; j < 4; ++j ) {
                    ComputeFloat tmp = crx_p * ( Sy0[j-1] + static_cast<ComputeFloat>( 0.5 ) * ( Sy1[j] - Sy0[j-1] )) ;
                    ComputeFloat tmp_reduction{};
                    const int iloc = ipo *  GPUClusterWithGCWidth + jpo + j;
                    tmp_reduction -= Sx1[0] * tmp;
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int i = 2; i < 5; ++i ) {
                        tmp_reduction -= ( Sx1[i - 1] - Sx0[i - 2] ) * tmp;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc + i * GPUClusterWithGCWidth], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }//*/


                // Jy

                for( unsigned int i = 0; i < 1; ++i ) {
                    const int    iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 , 2  ) + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; j++ ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int    iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 , 2  ) + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }//*/
                /*{
                    ComputeFloat tmp = cry_p * Sy1[0] * one_half;
                    const int    iloc = ( ipo ) * GPUClusterWithGCWidth + jpo;
                    ComputeFloat tmp_reduction{};
                    tmp_reduction -= Sx1[0] * tmp;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int i = 1; i < 5; ++i ) {
                        tmp_reduction -= tmp * ( Sx0[i-1] +  Sx1[i] );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + i * GPUClusterWithGCWidth], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }


                for( unsigned int j = 2; j < 5; ++j ) {
                     ComputeFloat tmp = cry_p * ( Sy1[j - 1] - Sy0[j - 2] ) * one_half;
                    const int    iloc = ( ipo ) * GPUClusterWithGCWidth + jpo + j;
                    ComputeFloat tmp_reduction{};
                    tmp_reduction -= Sx1[0] * tmp;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int i = 1; i < 5; ++i ) {
                        tmp_reduction -= tmp * ( Sx0[i-1] +  Sx1[i] );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + i * GPUClusterWithGCWidth], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }//*/


                // Jz

                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 , 2  ) + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * (  Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i]  ) +
                                                                                                                      Sy1[j] * (  Sx1[i] ) ) ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 , 2  ) + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] + Sx0[i] ) +
                                                                                                                      Sy1[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    }
                }//*/
                /*{
                    ComputeFloat tmp = crz_p * Sx1[0];
                    const int iloc = ipo * GPUClusterWithGCWidth + jpo;
                     atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( tmp *  Sy1[0]) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( tmp * (Sy0[j-1] * one_half + Sy1[j]) ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * GPUClusterWithGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( one_half * Sx0[i-1] + Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j-1] * ( one_half * Sx1[i] + Sx0[i-1] ) +
                                                                                                                      Sy1[j] * ( one_half * Sx0[i-1] + Sx1[i] ) ) ) );
                    }
                }//*/

            }

            __syncthreads();

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int local_x_scratch_space_coordinate = field_index / GPUClusterWithGCWidth;
                const unsigned int local_y_scratch_space_coordinate = field_index % GPUClusterWithGCWidth;

                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate;
                const unsigned int global_y_scratch_space_coordinate = global_y_scratch_space_coordinate_offset + local_y_scratch_space_coordinate;

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int global_memory_index = global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate;
                const unsigned int scratch_space_index = field_index; // local_x_scratch_space_coordinate * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + local_y_scratch_space_coordinate;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + /* We handle the FTDT/picsar */ not_spectral_ * global_x_scratch_space_coordinate], static_cast<double>( Jy_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[scratch_space_index] ) );
            }
        } // end DepositCurrent


        template <typename ComputeFloat,
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize>
        __global__ void
        // __launch_bounds__(kWorkgroupSize, 1)
        DepositCurrentAndDensity_2D_Order2( double *__restrict__ device_Jx,
                                            double *__restrict__ device_Jy,
                                            double *__restrict__ device_Jz,
                                            double *__restrict__ device_rho,
                                            int Jx_size,
                                            int Jy_size,
                                            int Jz_size,
                                            int rho_size,
                                            const double *__restrict__ device_particle_position_x,
                                            const double *__restrict__ device_particle_position_y,
                                            const double *__restrict__ device_particle_momentum_z,
                                            const short *__restrict__ device_particle_charge,
                                            const double *__restrict__ device_particle_weight,
                                            const int *__restrict__ device_bin_index,
                                            const double *__restrict__ device_invgf_,
                                            const int *__restrict__ device_iold_,
                                            const double *__restrict__ device_deltaold_,
                                            ComputeFloat inv_cell_volume,
                                            ComputeFloat dx_inv,
                                            ComputeFloat dy_inv,
                                            ComputeFloat dx_ov_dt,
                                            ComputeFloat dy_ov_dt,
                                            int          i_domain_begin,
                                            int          j_domain_begin,
                                            int          nprimy,
                                            int          not_spectral_ )
        {
            // TODO(Etienne M): refactor this function. Break it into smaller
            // pieces (lds init/store, coeff computation, deposition etc..)
            // TODO(Etienne M): __ldg could be used to slightly improve GDS load
            // speed. This would only have an effect on Nvidia cards as this
            // operation is a no op on AMD.
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x * gridDim.y;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y + y_cluster_coordinate; // The indexing order is: x * ywidth * zwidth + y * zwidth + z
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 2 /* 2D */ );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 2 /* 2D */ );

            // NOTE: We gain from the particles not being sorted inside a
            // cluster because it reduces the bank conflicts one gets when
            // multiple threads access the same part of the shared memory. Such
            // "conflicted" accesses are serialized !
            // NOTE: We use a bit to much LDS. For Jx, the first row could be
            // discarded, for Jy we could remove the first column.

            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 2 /* 2D */, 2 /* 2nd order interpolation */ );

            // NOTE: I tried having only one cache and reusing it. Doing that
            // requires you to iterate multiple time over the particle which is
            // possible but cost more bandwidth. The speedup was ~x0.92.
            __shared__ ReductionFloat Jx_scratch_space[kFieldScratchSpaceSize];
            __shared__ ReductionFloat Jy_scratch_space[kFieldScratchSpaceSize];
            __shared__ ReductionFloat Jz_scratch_space[kFieldScratchSpaceSize];
            __shared__ ReductionFloat rho_scratch_space[kFieldScratchSpaceSize];

            // Init the shared memory

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {
                Jx_scratch_space[field_index]  = static_cast<ReductionFloat>( 0.0 );
                Jy_scratch_space[field_index]  = static_cast<ReductionFloat>( 0.0 );
                Jz_scratch_space[field_index]  = static_cast<ReductionFloat>( 0.0 );
                rho_scratch_space[field_index] = static_cast<ReductionFloat>( 0.0 );
            }

            __syncthreads();

            const unsigned int particle_count = device_bin_index[bin_count - 1];

            // This workgroup has to process distance(last_particle,
            // first_particle) particles
            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 :
                                                                                     device_bin_index[workgroup_dedicated_bin_index - 1];
            const unsigned int last_particle  = device_bin_index[workgroup_dedicated_bin_index];

            for( unsigned int particle_index = first_particle + thread_index_offset;
                 particle_index < last_particle;
                 particle_index += loop_stride ) {
                const ComputeFloat invgf                  = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__ iold        = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];
                // double DSx[5];
                // double DSy[5];

                // Variable declaration & initialization
                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // Locate the particle on the primal grid at former time-step & calculate coeff. S0
                {
                    const ComputeFloat delta  = deltaold[0 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sx0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[4] = static_cast<ComputeFloat>( 0.0 );
                }
                {
                    const ComputeFloat delta  = deltaold[1 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sy0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[4] = static_cast<ComputeFloat>( 0.0 );
                }

                // Locate the particle on the primal grid at current time-step & calculate coeff. S1
                {
                    // const int    ip             = static_cast<int>( xpn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat xpn      = static_cast<ComputeFloat>( device_particle_position_x[particle_index] ) * dx_inv;
                    const int          ip       = std::round( xpn );
                    const int          ipo      = iold[0 * particle_count];
                    const int          ip_m_ipo = ip - ipo - i_domain_begin;
                    const ComputeFloat delta    = xpn - static_cast<ComputeFloat>( ip );
                    const ComputeFloat delta2   = delta * delta;

                    Sx1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sx1[2] = 0.0; // Always set below
                    Sx1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sx1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }
                {
                    // const int    jp             = static_cast<int>( ypn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat ypn      = static_cast<ComputeFloat>( device_particle_position_y[particle_index] ) * dy_inv;
                    const int          jp       = std::round( ypn );
                    const int          jpo      = iold[1 * particle_count];
                    const int          jp_m_jpo = jp - jpo - j_domain_begin;
                    const ComputeFloat delta    = ypn - static_cast<ComputeFloat>( jp );
                    const ComputeFloat delta2   = delta * delta;

                    Sy1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sy1[2] = 0.0; // Always set below
                    Sy1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sy1[jp_m_jpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy1[jp_m_jpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy1[jp_m_jpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }

                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( 1.0 / 3.0 ) * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;

                // Jx

                ComputeFloat tmpJx[5]{};

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( static_cast<ComputeFloat>( 0.5 ) * ( Sy1[0] - Sy0[0] ) );
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmpJx[0] ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + static_cast<ComputeFloat>( 0.5 ) * ( Sy1[j] - Sy0[j] ) );
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc + j], static_cast<ReductionFloat>( tmpJx[j] ) );
                    }
                }

                // Jy

                for( unsigned int i = 0; i < 1; ++i ) {
                    const int    iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; j++ ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int    iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }

                // Jz

                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] /* + Sx0[i] */ ) +
                                                                                                                      Sy1[j] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) ) ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] + Sx0[i] ) +
                                                                                                                      Sy1[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    }
                }

                // Rho
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    atomic::LDS::AddNoReturn( &rho_scratch_space[iloc], static_cast<ReductionFloat>( charge_weight * ( Sx1[0] * Sy1[i] ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &rho_scratch_space[iloc + j], static_cast<ReductionFloat>( charge_weight * ( Sx1[0] * Sy1[j] ) ) );
                    }
                }

                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + jpo;
                    atomic::LDS::AddNoReturn( &rho_scratch_space[iloc], static_cast<ReductionFloat>( charge_weight * ( Sx1[i] * Sy1[0] ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &rho_scratch_space[iloc + j], static_cast<ReductionFloat>( charge_weight * ( Sx1[i] * Sy1[j] ) ) );
                    }
                }
            }

            __syncthreads();

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int local_x_scratch_space_coordinate = field_index / Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ );
                const unsigned int local_y_scratch_space_coordinate = field_index % Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ );

                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate;
                const unsigned int global_y_scratch_space_coordinate = global_y_scratch_space_coordinate_offset + local_y_scratch_space_coordinate;

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int global_memory_index = global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate;
                const unsigned int scratch_space_index = field_index; // local_x_scratch_space_coordinate * Params::getGPUClusterWithGhostCellWidth( 2 /* 2D */, 2 /* 2nd order interpolation */ ) + local_y_scratch_space_coordinate;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + /* We handle the FTDT/picsar */ not_spectral_ * global_x_scratch_space_coordinate], static_cast<double>( Jy_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_rho[global_memory_index], static_cast<double>( rho_scratch_space[scratch_space_index] ) );
            }
        }
    } // namespace kernel


    //static inline
    void
    currentDepositionKernel2D( double *__restrict__ host_Jx,
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
        SMILEI_ASSERT( Params::getGPUClusterWidth( 2 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        // NOTE:
        // This cluster is very strongly bound by atomic operations in LDS (shared memory)
        // TODO(Etienne M): Find a way to lessen the atomic usage

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), 1 };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // NOTE: On cards lacking hardware backed Binary64 atomic operations,
        // falling back to Binary32 (supposing hardware support for atomic
        // operations) can lead to drastic performance improvement.
        // One just need to assign 'float' to ReductionFloat.
        //
        using ComputeFloat   = double;
        using ReductionFloat = double;

	auto KernelFunction = kernel::DepositCurrentDensity_2D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;
#if defined ( __HIP__ ) 
        hipLaunchKernelGGL( KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, // Shared memory
                            0, // Stream
                            // Kernel arguments
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jx ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jy ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jz ),
                            Jx_size, Jy_size, Jz_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv, dy_inv,
                            dx_ov_dt, dy_ov_dt,
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_ );

        checkHIPErrors( ::hipDeviceSynchronize() );
#elif defined (  __NVCC__ )
	KernelFunction <<<
                            kGridDimension,
                            kBlockDimension,
                            0, // Shared memory
                            0 // Stream
                       >>>
                       (
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jx ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jy ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jz ),
                            Jx_size, Jy_size, Jz_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv, dy_inv,
                            dx_ov_dt, dy_ov_dt,
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif
    }

    //static inline 
    void
    currentAndDensityDepositionKernel2D( double *__restrict__ host_Jx,
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
        SMILEI_ASSERT( Params::getGPUClusterWidth( 2 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        // NOTE:
        // This cluster is very strongly bound by atomic operations in LDS (shared memory)
        // TODO(Etienne M): Find a way to lessen the atomic usage

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), 1 };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // NOTE: On cards lacking hardware backed Binary64 atomic operations,
        // falling back to Binary32 (supposing hardware support for atomic
        // operations) can lead to drastic performance improvement.
        // One just need to assign 'float' to ReductionFloat.
        //
        using ComputeFloat   = double;
        using ReductionFloat = double;
        auto KernelFunction = kernel::DepositCurrentAndDensity_2D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;
#if defined ( __HIP__ ) 
        hipLaunchKernelGGL( KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, // Shared memory
                            0, // Stream
                            // Kernel arguments
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jx ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jy ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jz ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_rho ),
                            Jx_size, Jy_size, Jz_size, rho_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv, dy_inv,
                            dx_ov_dt, dy_ov_dt,
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_ );

        checkHIPErrors( ::hipDeviceSynchronize() );
#elif defined (  __NVCC__ )
        KernelFunction <<<                                                                                                             
                            kGridDimension,                                                                                            
                            kBlockDimension,                                                                                           
                            0, // Shared memory                                                                                        
                            0 // Stream                                                                                                
                       >>>                                                                                                             
                       (                                                                                                               
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jx ),                               
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jy ),                               
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jz ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_rho ),
                            Jx_size, Jy_size, Jz_size, rho_size,
                            device_particle_position_x,                                                                                
                            device_particle_position_y,                                                                                
                            device_particle_momentum_z,                                                                                
                            device_particle_charge,                                                                                    
                            device_particle_weight,                                                                                    
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),                        
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),                           
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),                            
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),                        
                            inv_cell_volume,                                                                                           
                            dx_inv, dy_inv,                                                                                            
                            dx_ov_dt, dy_ov_dt,                                                                                        
                            i_domain_begin, j_domain_begin,                                                                            
                            nprimy,                                                                                                    
                            not_spectral_                                                                                               
                       );                                                                                                              
        checkHIPErrors( ::cudaDeviceSynchronize() );  
#endif 
    }

} // namespace hipcuda

//#endif

//! Project global current densities (EMfields->Jx_/Jy_/Jz_)
//!
//extern "C" void
//currentDepositionKernel2D( double *__restrict__ host_Jx,
//                         double *__restrict__ host_Jy,
//                         double *__restrict__ host_Jz,
//                         int Jx_size,
//                         int Jy_size,
//                         int Jz_size,
//                         const double *__restrict__ device_particle_position_x,
//                         const double *__restrict__ device_particle_position_y,
//                         const double *__restrict__ device_particle_momentum_z,
//                         const short *__restrict__ device_particle_charge,
//                         const double *__restrict__ device_particle_weight,
//                         const int *__restrict__ host_bin_index,
//                         unsigned int x_dimension_bin_count,
//                         unsigned int y_dimension_bin_count,
//                         const double *__restrict__ host_invgf_,
//                         const int *__restrict__ host_iold_,
//                         const double *__restrict__ host_deltaold_,
//                         double inv_cell_volume,
//                         double dx_inv,
//                         double dy_inv,
//                         double dx_ov_dt,
//                         double dy_ov_dt,
//                         int    i_domain_begin,
//                         int    j_domain_begin,
//                         int    nprimy,
//                         int    not_spectral_ )
//{
//    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
//    naive:: // the naive, OMP version serves as a reference along with the CPU version
//    #else
//    cudahip::
//    #endif
//        currentDepositionKernel2D( host_Jx, host_Jy, host_Jz,
//                                 Jx_size, Jy_size, Jz_size,
//                                 device_particle_position_x, device_particle_position_y,
//                                 device_particle_momentum_z,
//                                 device_particle_charge,
//                                 device_particle_weight,
//                                 host_bin_index,
//                                 x_dimension_bin_count,
//                                 y_dimension_bin_count,
//                                 host_invgf_,
//                                 host_iold_, host_deltaold_,
//                                 inv_cell_volume,
//                                 dx_inv, dy_inv,
//                                 dx_ov_dt, dy_ov_dt,
//                                 i_domain_begin, j_domain_begin,
//                                 nprimy,
//                                 not_spectral_ );
//}
//
////! Project global current and charge densities (EMfields->Jx_/Jy_/Jz_/rho_)
////!
//extern "C" void
//currentAndDensityDepositionKernel( double *__restrict__ host_Jx,
//                                   double *__restrict__ host_Jy,
//                                   double *__restrict__ host_Jz,
//                                   double *__restrict__ host_rho,
//                                   int Jx_size,
//                                   int Jy_size,
//                                   int Jz_size,
//                                   int rho_size,
//                                   const double *__restrict__ device_particle_position_x,
//                                   const double *__restrict__ device_particle_position_y,
//                                   const double *__restrict__ device_particle_momentum_z,
//                                   const short *__restrict__ device_particle_charge,
//                                   const double *__restrict__ device_particle_weight,
//                                   const int *__restrict__ host_bin_index,
//                                   unsigned int x_dimension_bin_count,
//                                   unsigned int y_dimension_bin_count,
//                                   const double *__restrict__ host_invgf_,
//                                   const int *__restrict__ host_iold_,
//                                   const double *__restrict__ host_deltaold_,
//                                   double inv_cell_volume,
//                                   double dx_inv,
//                                   double dy_inv,
//                                   double dx_ov_dt,
//                                   double dy_ov_dt,
//                                   int    i_domain_begin,
//                                   int    j_domain_begin,
//                                   int    nprimy,
//                                   int    not_spectral_ )
//{
//    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
//    naive:: // the naive, OMP version serves as a reference along with the CPU version
//    #else
//    cudahip::
//    #endif
//        currentAndDensityDepositionKernel( host_Jx, host_Jy, host_Jz, host_rho,
//                                           Jx_size, Jy_size, Jz_size, rho_size,
//                                           device_particle_position_x, device_particle_position_y,
//                                           device_particle_momentum_z,
//                                           device_particle_charge,
//                                           device_particle_weight,
//                                           host_bin_index,
//                                           x_dimension_bin_count,
//                                           y_dimension_bin_count,
//                                           host_invgf_,
//                                           host_iold_, host_deltaold_,
//                                           inv_cell_volume,
//                                           dx_inv, dy_inv,
//                                           dx_ov_dt, dy_ov_dt,
//                                           i_domain_begin, j_domain_begin,
//                                           nprimy,
//                                           not_spectral_ );
//}

