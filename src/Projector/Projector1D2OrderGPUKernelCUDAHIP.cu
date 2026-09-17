

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

namespace cudahip1d {
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
        #if defined( __gfx90a__ ) ||  defined (__gfx942__)
                    ::unsafeAtomicAdd( a_pointer, a_value );
        #else
                    ::atomicAdd( a_pointer, a_value );
        #endif
                }

                __device__ void
                AddNoReturn( double *a_pointer, double a_value )
                {
        #if defined( __gfx90a__ ) ||  defined (__gfx942__)
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
        #if defined( __gfx90a__ ) ||  defined (__gfx942__)
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
        DepositCurrentDensity_1D_Order2( double *__restrict__ device_Jx,
                                         double *__restrict__ device_Jy,
                                         double *__restrict__ device_Jz,
                                         int Jx_size,
                                         int Jy_size,
                                         int Jz_size,
                                         const double *__restrict__ device_particle_position_x,
                                         const double *__restrict__ device_particle_momentum_y,
                                         const double *__restrict__ device_particle_momentum_z,
                                         const short *__restrict__ device_particle_charge,
                                         const double *__restrict__ device_particle_weight,
                                         const int *__restrict__ device_bin_index,
                                         const double *__restrict__ device_invgf_,
                                         const int *__restrict__ device_iold_,
                                         const double *__restrict__ device_deltaold_,
                                         ComputeFloat inv_cell_volume,
                                         ComputeFloat dx_inv,
                                         ComputeFloat dx_ov_dt,
                                         int          i_domain_begin,
                                         int          not_spectral_,
                                         unsigned int oversize_,
                                         bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate;
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 1 /* 1D */ );
            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 1 /* 1D */, 2 /* 2nd order interpolation */ );

            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 1 /* 1D */, 2 /* 2nd order interpolation */ );

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

            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 : device_bin_index[workgroup_dedicated_bin_index - 1];
            const unsigned int last_particle  = device_bin_index[workgroup_dedicated_bin_index];

            // The loop order is different depending on cell sorting
            unsigned int stride, start_thread, stop_thread;
            if( cell_sorting ) {
                // With cell sorting, each thread should process close-by particles
                // to reduce atomics. This uses more cache, but is still better
                const unsigned int npart_thread = last_particle > first_particle ? ( last_particle - first_particle - 1 ) / workgroup_size + 1 : 0;
                start_thread = first_particle + threadIdx.x * npart_thread;
                stop_thread = std::min( { start_thread + npart_thread, last_particle } );
                stride  = 1;
            } else {
                // Without cell sorting, we keep the standard loops as particles
                // are not ordered so that atomics are naturally rare
                start_thread = first_particle + threadIdx.x;
                stop_thread = last_particle;
                stride = workgroup_size;
            }
            for( unsigned int particle_index = start_thread; particle_index < stop_thread; particle_index += stride ) {
                const ComputeFloat invgf                  = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__ iold        = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * static_cast<ComputeFloat>( device_particle_momentum_y[particle_index] ) * invgf;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                ComputeFloat Sx0[3];
                ComputeFloat Sx1[5];

                // Variable declaration & initialization
                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // Locate the particle on the primal grid at former time-step & calculate coeff. S0
                {
                    const ComputeFloat delta  = deltaold[0 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sx0[0] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[1] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[2] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }
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

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);

                // Jx
                ComputeFloat tmpJx = 0.0; 
                // i=1
                tmpJx += crx_p * ( - Sx1[0] ); 
                atomic::LDS::AddNoReturn( &Jx_scratch_space[ipo + 1], static_cast<ReductionFloat>( tmpJx ) );
                for( unsigned int i = 2; i < 5; ++i ) {
                    const int    iloc = i + ipo;
                    tmpJx = tmpJx + crx_p * (Sx0[i-2] - Sx1[i-1]); 
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmpJx ) );
                }
                // Jy & Jz
                //i=0
                {
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[ipo], static_cast<ReductionFloat>( cry_p * 0.5 * Sx1[0] ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[ipo], static_cast<ReductionFloat>( crz_p * 0.5 * Sx1[0] ) );
                }
                for( unsigned int i = 1; i < 4; ++i ) {
                    const int    iloc = i + ipo;
                    double temp = 0.5 * (Sx0[i-1] + Sx1[i]);
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( cry_p * temp ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * temp ) );
                }
                //i=4
                {
                    const int    iloc = 4 + ipo;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( cry_p * 0.5 * Sx1[4] ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * 0.5 * Sx1[4] ) );
                }
            } // particle_index

            __syncthreads();

            for( unsigned int field_index = thread_index_offset; field_index < kFieldScratchSpaceSize; field_index += workgroup_size ) {
                const unsigned int local_x_scratch_space_coordinate = field_index % GPUClusterWithGCWidth; // /GPUClusterWithGCWidth
                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate+(oversize_-2);

                const unsigned int global_memory_index = global_x_scratch_space_coordinate;
                const unsigned int scratch_space_index = field_index; // local_x_scratch_space_coordinate * GPUClusterWithGCWidth + local_y_scratch_space_coordinate;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index /*+  not_spectral * global_x_scratch_space_coordinate*/], static_cast<double>( Jy_scratch_space[scratch_space_index] ) ); //  We handle the FTDT/picsar 
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[scratch_space_index] ) );
            }
        } // end DepositCurrent


        template <typename ComputeFloat,
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize>
        __global__ void
        // __launch_bounds__(kWorkgroupSize, 1)
        DepositCurrentAndDensity_1D_Order2( double *__restrict__ device_Jx,
                                            double *__restrict__ device_Jy,
                                            double *__restrict__ device_Jz,
                                            double *__restrict__ device_rho,
                                            int Jx_size,
                                            int Jy_size,
                                            int Jz_size,
                                            int rho_size,
                                            const double *__restrict__ device_particle_position_x,
                                            const double *__restrict__ device_particle_momentum_y,
                                            const double *__restrict__ device_particle_momentum_z,
                                            const short *__restrict__ device_particle_charge,
                                            const double *__restrict__ device_particle_weight,
                                            const int *__restrict__ device_bin_index,
                                            const double *__restrict__ device_invgf_,
                                            const int *__restrict__ device_iold_,
                                            const double *__restrict__ device_deltaold_,
                                            ComputeFloat inv_cell_volume,
                                            ComputeFloat dx_inv,
                                            ComputeFloat dx_ov_dt,
                                            int          i_domain_begin,
                                            int          not_spectral_,
                                            unsigned int oversize_,
                                            bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate ; 
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 1 /* 1D */ );

            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 1 /* 1D */, 2 /* 2nd order interpolation */ );
            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 1 /* 1D */, 2 /* 2nd order interpolation */ );

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

            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 : device_bin_index[workgroup_dedicated_bin_index - 1];
            const unsigned int last_particle  = device_bin_index[workgroup_dedicated_bin_index];

            unsigned int stride, start_thread, stop_thread;
            if( cell_sorting ) {
                const unsigned int npart_thread = last_particle > first_particle ? ( last_particle - first_particle - 1 ) / workgroup_size + 1 : 0;
                start_thread = first_particle + threadIdx.x * npart_thread;
                stop_thread = std::min( { start_thread + npart_thread, last_particle } );
                stride  = 1;
            } else {
                start_thread = first_particle + threadIdx.x;
                stop_thread = last_particle;
                stride = workgroup_size;
            }

            for( unsigned int particle_index = start_thread; particle_index < stop_thread; particle_index += stride ) {
                const ComputeFloat                  invgf = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__        iold = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * static_cast<ComputeFloat>( device_particle_momentum_y[particle_index] ) * invgf;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                ComputeFloat Sx0[3];
                ComputeFloat Sx1[5];

                // Variable declaration & initialization
                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // Locate the particle on the primal grid at former time-step & calculate coeff. S0
                {
                    const ComputeFloat delta  = deltaold[0 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sx0[0] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[1] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[2] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
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

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);

                // Jx
                ComputeFloat tmpJx = 0.0; 
                // i=1
                tmpJx += crx_p * ( - Sx1[0] ); 
                atomic::LDS::AddNoReturn( &Jx_scratch_space[ipo + 1], static_cast<ReductionFloat>( tmpJx ) );
                for( unsigned int i = 2; i < 5; ++i ) {
                    const int    iloc = i + ipo;
                    tmpJx = tmpJx + crx_p * (Sx0[i-2] - Sx1[i-1]); 
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[iloc], static_cast<ReductionFloat>( tmpJx ) );
                }
                // Jy & Jz
                //i=0
                {
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[ipo], static_cast<ReductionFloat>( cry_p * 0.5 * Sx1[0] ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[ipo], static_cast<ReductionFloat>( crz_p * 0.5 * Sx1[0] ) );
                }
                for( unsigned int i = 1; i < 4; ++i ) {
                    const int    iloc = i + ipo;
                    double temp = 0.5 * (Sx0[i-1] + Sx1[i]);
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( cry_p * temp ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * temp ) );
                }
                //i=4
                {
                    const int    iloc = 4 + ipo;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[iloc], static_cast<ReductionFloat>( cry_p * 0.5 * Sx1[4] ) );
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[iloc], static_cast<ReductionFloat>( crz_p * 0.5 * Sx1[4] ) );
                }

                // Rho
                for( unsigned int i = 0; i < 5; ++i ) {
                    const int iloc = i + ipo;
                    atomic::LDS::AddNoReturn( &rho_scratch_space[iloc], static_cast<ReductionFloat>( charge_weight * Sx1[i] ) );
                }

            } // particle_index

            __syncthreads();

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {

                const unsigned int local_x_scratch_space_coordinate = field_index % GPUClusterWithGCWidth;
                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate+(oversize_-2);

                const unsigned int global_memory_index = global_x_scratch_space_coordinate;
                const unsigned int scratch_space_index = field_index;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[scratch_space_index] ) );
                //atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + /* We handle the FTDT/picsar */ not_spectral * global_x_scratch_space_coordinate], static_cast<double>( Jy_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index], static_cast<double>( Jy_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[scratch_space_index] ) );
                atomic::GDS::AddNoReturn( &device_rho[global_memory_index], static_cast<double>( rho_scratch_space[scratch_space_index] ) );
            }
        }
    } // namespace kernel


    //static inline
    void
    currentDepositionKernel1D( double *__restrict__ host_Jx,
                             double *__restrict__ host_Jy,
                             double *__restrict__ host_Jz,
                             int Jx_size,
                             int Jy_size,
                             int Jz_size,
                             const double *__restrict__ device_particle_position_x,
                             const double *__restrict__ device_particle_momentum_y,
                             const double *__restrict__ device_particle_momentum_z,
                             const short *__restrict__ device_particle_charge,
                             const double *__restrict__ device_particle_weight,
                             const int *__restrict__ host_bin_index,
                             unsigned int x_dimension_bin_count,
                             const double *__restrict__ host_invgf_,
                             const int *__restrict__ host_iold_,
                             const double *__restrict__ host_deltaold_,
                             double inv_cell_volume,
                             double dx_inv,
                             double dx_ov_dt,
                             int    i_domain_begin,
                             int    not_spectral_,
                             unsigned int oversize_,
                             bool   cell_sorting )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 1 /* 1D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        // NOTE: This cluster is very strongly bound by atomic operations in LDS (shared memory)

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), 1, 1 };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // NOTE: On cards lacking hardware backed Binary64 atomic operations,
        // falling back to Binary32 (supposing hardware support for atomic
        // operations) can lead to drastic performance improvement.
        // One just need to assign 'float' to ReductionFloat.
        //
        using ComputeFloat   = double;
        using ReductionFloat = double;

	auto KernelFunction = kernel::DepositCurrentDensity_1D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;
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
                            device_particle_momentum_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv,
                            dx_ov_dt,
                            i_domain_begin,
                            not_spectral_,
                            oversize_,
                            cell_sorting );

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
                            device_particle_momentum_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv,
                            dx_ov_dt,
                            i_domain_begin,
                            not_spectral_,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif
    }

    //static inline 
    void
    currentAndDensityDepositionKernel1D( double *__restrict__ host_Jx,
                                       double *__restrict__ host_Jy,
                                       double *__restrict__ host_Jz,
                                       double *__restrict__ host_rho,
                                       int Jx_size,
                                       int Jy_size,
                                       int Jz_size,
                                       int rho_size,
                                       const double *__restrict__ device_particle_position_x,
                                       const double *__restrict__ device_particle_momentum_y,
                                       const double *__restrict__ device_particle_momentum_z,
                                       const short *__restrict__ device_particle_charge,
                                       const double *__restrict__ device_particle_weight,
                                       const int *__restrict__ host_bin_index,
                                       unsigned int x_dimension_bin_count,
                                       const double *__restrict__ host_invgf_,
                                       const int *__restrict__ host_iold_,
                                       const double *__restrict__ host_deltaold_,
                                       double inv_cell_volume,
                                       double dx_inv,
                                       double dx_ov_dt,
                                       int    i_domain_begin,
                                       int    not_spectral_,
                                       unsigned int oversize_,
                                       bool cell_sorting )
    {
        // & because one  1D ; 2 because of 2nd order interpolation
        SMILEI_ASSERT( Params::getGPUClusterWidth( 1 ) != -1 && 
                       Params::getGPUClusterGhostCellBorderWidth( 2 ) != -1 );

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), 1, 1 };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        using ComputeFloat   = double;
        using ReductionFloat = double;
        auto KernelFunction = kernel::DepositCurrentAndDensity_1D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;
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
                            device_particle_momentum_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv,
                            dx_ov_dt,
                            i_domain_begin,
                            not_spectral_,
                            oversize_,
                            cell_sorting );

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
                            device_particle_momentum_y,
                            device_particle_momentum_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv,
                            dx_ov_dt,
                            i_domain_begin,
                            not_spectral_,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif 
    }

} // namespace cudahip1D


