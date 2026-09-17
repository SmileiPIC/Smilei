//! HIP CUDA implementation
//!
//! ====================================================================
//! Projector3D2OrderGPUKernelCUDAHIP 
//! ====================================================================
//! Increase Projector performance thanks to Nsight-Compute report
//! --------------------------------------------------------------------
//! Mixed precision
//! --------------------------------------------------------------------
//! To increase the computational throughput some calculation can use
//! float instead of double without degrading numerical heating.
//! For example using 'float' for reduction allow speed-up which can reduce
//! contention. In top of that :
//! Adds a ComputePositionFloat type (typically double) used
//! exclusively for the position → delta calculation, which involves
//! catastrophic cancellation (xpn - round(xpn)). All other
//! computations (shape functions, weights, deposition) use
//! ComputeFloat (can be float for register savings).
//!
//! Precision map:
//!   ComputePositionFloat (double)  — position × dx_inv, xpn - ip
//!   ComputeFloat         (float)   — shape functions, weights, J deposition
//!   ReductionFloat       (float)   — shared memory atomics
//!   double                         — global memory flush (always)
//! --------------------------------------------------------------------
//! TUNING GUIDE
//! --------------------------------------------------------------------
//! In order to maximize the occupancy and the correct use of GPU, some
//! parameters have been add to tune the projector in accordance with the GPU
//! capability.
//!
//! Goals :
//!   - Reduce contention
//!   - Increase or deacrese the use of LDS/Shared Memory to adjust the number of
//!     Blocks per Stream Multiprocessor
//!   - Increase or deacrese the use of register to adjust the number of
//!     Blocks per Stream Multiprocessor
//!
//!   kWorkgroupSize   — Threads per block (typically 128 in 3D)
//!   kNBufferCopies   — Replicated shared memory buffers (1 = off) to reduce contention (useful in 2D !)
//!   kShmemPad        — Shared memory stride padding (0 = off) to reduce bank conflict
//!   kMinBlocksPerSM  — Minimum blocks per SM for __launch_bounds__ to force a number of block per SM
//!                      to increase occupancy if LDS or Register allow it
//!
//! Bonuses with mixed precision :
//!
//! Type aliases:
//!   ComputeFloat         = float  → saves ~27 registers (double→float on shapes)
//!   ComputePositionFloat = double → preserves position precision
//!   ReductionFloat       = float  → shared memory atomics (unchanged)
//!
//! Tested configurations on GH200 (228 KB shmem/SM, 65536 regs/SM):
//!
//!   Config A: copies=1, pad=0, minBlocks=8  → 89 ms
//!   Config B: copies=2, pad=0, minBlocks=8  → 89 ms
//!   Config C: copies=1, pad=2, minBlocks=8  → 85 ms (if bank conflicts dominate)
//!   Config D: copies=2, pad=2, minBlocks=6  → 101 ms (shmem limited)
//!
//! Bank conflict padding math (for GCWidth=9):
//!   pad=0 → width= 9, strides (1, 9, 81),  81 mod 32 = 17 (ok-ish)
//!   pad=1 → width=10, strides (1,10,100), 100 mod 32 =  4 (BAD)
//!   pad=2 → width=11, strides (1,11,121), 121 mod 32 = 25 (GOOD)
//!   pad=3 → width=12, strides (1,12,144), 144 mod 32 = 16 (BAD)
//! 
//! All theses configuration have been tested and stress with Nsight-Compute
//! ====================================================================

// TODO(Etienne M): The makefile does not recognise this file and doesn't compute
// it's dependencies. If you make a modification in one of the header this file
// includes, you must `touch` this file. IF you dont do that you'll have ABI/ODR
// issues (!).

#if defined( SMILEI_ACCELERATOR_GPU )


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


        // ====================================================================
        // Current deposition kernel (Jx, Jy, Jz) — without rho
        // ====================================================================

        template <typename ComputeFloat,
                  typename ComputePositionFloat,       // High-precision type for positions
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize,
                  std::size_t kNBufferCopies,
                  std::size_t kShmemPad,
                  std::size_t kMinBlocksPerSM>
        __global__ void
        __launch_bounds__(kWorkgroupSize, kMinBlocksPerSM)
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
                                         ComputePositionFloat dx_inv,  // Position precision
                                         ComputePositionFloat dy_inv,  // Position precision
                                         ComputeFloat dx_ov_dt,
                                         ComputeFloat dy_ov_dt,
                                         int          i_domain_begin,
                                         int          j_domain_begin,
                                         int          nprimy,
                                         int          not_spectral_,
                                         unsigned int oversize_,
                                         bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize;
            const unsigned int bin_count      = gridDim.x * gridDim.y;

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y + y_cluster_coordinate;
            const unsigned int thread_index_offset           = threadIdx.x;

            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 2 );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 2 );

            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 2, 2 );
            const int PaddedGCWidth = GPUClusterWithGCWidth + static_cast<int>( kShmemPad ); // To reduce Bank Conflict

            static constexpr unsigned int kLogicalGCWidth   = Params::getGPUClusterWithGhostCellWidth( 2, 2 );
            static constexpr unsigned int kPaddedGCWidth    = kLogicalGCWidth + kShmemPad;
            static constexpr unsigned int kPaddedFieldSize  = kPaddedGCWidth * kPaddedGCWidth;

            // NOTE: I tried having only one cache and reusing it. Doing that
            // requires you to iterate multiple time over the particle which is
            // possible but cost more bandwidth. The speedup was ~x0.92.

            __shared__ ReductionFloat Jx_scratch_space[kNBufferCopies][kPaddedFieldSize]; // Copie to reduce contention
            __shared__ ReductionFloat Jy_scratch_space[kNBufferCopies][kPaddedFieldSize];
            __shared__ ReductionFloat Jz_scratch_space[kNBufferCopies][kPaddedFieldSize];
            //extern __shared__ ReductionFloat J_scratch_space[];
            //ReductionFloat *Jx_scratch_space = J_scratch_space;
            //ReductionFloat *Jy_scratch_space = (ReductionFloat*)&Jx_scratch_space[kFieldScratchSpaceSize];// ReductionFloat *Jy_scratch_space = &J_scratch_space[1*kFieldScratchSpaceSize];
            //ReductionFloat *Jz_scratch_space = (ReductionFloat*)&Jy_scratch_space[kFieldScratchSpaceSize];// ReductionFloat *Jz_scratch_space = &J_scratch_space[2*kFieldScratchSpaceSize];

            const unsigned int buffer_copy_id = threadIdx.x % kNBufferCopies;

            for( unsigned int copy = 0; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
                    Jx_scratch_space[copy][field_index] = static_cast<ReductionFloat>( 0.0 );
                    Jy_scratch_space[copy][field_index] = static_cast<ReductionFloat>( 0.0 );
                    Jz_scratch_space[copy][field_index] = static_cast<ReductionFloat>( 0.0 );
                }
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
                const ComputeFloat invgf                  = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__ iold        = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];

                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // ── S0 coefficients ──────────────────────────────────
                // deltaold values are already deltas in [-0.5, 0.5],
                // stored as double in memory. Safe to cast to ComputeFloat
                // (no catastrophic cancellation here).
                {
                    const ComputeFloat delta  = static_cast<ComputeFloat>( deltaold[0 * particle_count] );  // Explicit cast
                    const ComputeFloat delta2 = delta * delta;
                    Sx0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[4] = static_cast<ComputeFloat>( 0.0 );
                }
                {
                    const ComputeFloat delta  = static_cast<ComputeFloat>( deltaold[1 * particle_count] );  // Explicit cast
                    const ComputeFloat delta2 = delta * delta;
                    Sy0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[4] = static_cast<ComputeFloat>( 0.0 );
                }

                // ── S1 coefficients ──────────────────────────────────
                // CRITICAL: position × dx_inv and xpn - ip MUST use
                // ComputePositionFloat (double) to avoid catastrophic
                // cancellation. The result (delta) is in [-0.5, 0.5]
                // and safely cast to ComputeFloat for shape functions.
                {
                    const ComputePositionFloat xpn   = static_cast<ComputePositionFloat>( device_particle_position_x[particle_index] ) * dx_inv;  // high precision
                    const int          ip            = std::round( xpn );                                                                          // round in high precision
                    const int          ipo           = iold[0 * particle_count];
                    const int          ip_m_ipo      = ip - ipo - i_domain_begin;
                    const ComputeFloat delta         = static_cast<ComputeFloat>( xpn - static_cast<ComputePositionFloat>( ip ) );                 // subtract in high precision, then narrow
                    const ComputeFloat delta2        = delta * delta;

                    Sx1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[1] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sx1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }
                {
                    const ComputePositionFloat ypn   = static_cast<ComputePositionFloat>( device_particle_position_y[particle_index] ) * dy_inv;  // high precision
                    const int          jp            = std::round( ypn );                                                                          // round in high precision
                    const int          jpo           = iold[1 * particle_count];
                    const int          jp_m_jpo      = jp - jpo - j_domain_begin;
                    const ComputeFloat delta         = static_cast<ComputeFloat>( ypn - static_cast<ComputePositionFloat>( jp ) );                 // subtract in high precision, then narrow
                    const ComputeFloat delta2        = delta * delta;

                    Sy1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[1] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sy1[jp_m_jpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy1[jp_m_jpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy1[jp_m_jpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }

                // ── Deposition weights ───────────────────────────────
                // All in ComputeFloat (no cancellation risk)
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( 1.0 / 3.0 ) * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;

                // ── Jx deposition ────────────────────────────────────
                ComputeFloat tmpJx[5]{};
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( static_cast<ComputeFloat>( 0.5 ) * ( Sy1[0] - Sy0[0] ) );
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( tmpJx[0] ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + static_cast<ComputeFloat>( 0.5 ) * ( Sy1[j] - Sy0[j] ) );
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmpJx[j] ) );
                    }
                }

                // ── Jy deposition ────────────────────────────────────
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int    iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; j++ ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int    iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }

                // ── Jz deposition ────────────────────────────────────
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] ) +
                                                                                                                      Sy1[j] * ( Sx1[i] ) ) ) );
                    }
                }
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] + Sx0[i] ) +
                                                                                                                      Sy1[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    }
                }
            }

            __syncthreads();

            // Reduction (no-op when kNBufferCopies=1)
            for( unsigned int copy = 1; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
                    Jx_scratch_space[0][field_index] += Jx_scratch_space[copy][field_index];
                    Jy_scratch_space[0][field_index] += Jy_scratch_space[copy][field_index];
                    Jz_scratch_space[0][field_index] += Jz_scratch_space[copy][field_index];
                }
            }

            if( kNBufferCopies > 1 ) {
                __syncthreads();
            }

            // Flush shared memory to global memory
            for( unsigned int field_index = thread_index_offset;
                 field_index < kPaddedFieldSize;
                 field_index += workgroup_size ) {

                const unsigned int local_x = field_index / PaddedGCWidth;
                const unsigned int local_y = field_index % PaddedGCWidth;

                if( kShmemPad > 0 ) {
                    if( local_x >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_y >= static_cast<unsigned int>(GPUClusterWithGCWidth) ) {
                        continue;
                    }
                }

                const unsigned int global_x = global_x_scratch_space_coordinate_offset + local_x + (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;
                const unsigned int global_y = global_y_scratch_space_coordinate_offset + local_y + (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;

                const unsigned int global_memory_index = global_x * nprimy + global_y;

                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + not_spectral_ * global_x], static_cast<double>( Jy_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[0][field_index] ) );
            }
        } // end DepositCurrentDensity_2D_Order2

        // ====================================================================
        // Current + Density deposition kernel (Jx, Jy, Jz, rho)
        // ====================================================================

        template <typename ComputeFloat,
                  typename ComputePositionFloat,       // NEW
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize,
                  std::size_t kNBufferCopies,
                  std::size_t kShmemPad,
                  std::size_t kMinBlocksPerSM>
        __global__ void
        __launch_bounds__(kWorkgroupSize, kMinBlocksPerSM)
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
                                            ComputePositionFloat dx_inv,  // Position precision
                                            ComputePositionFloat dy_inv,  // Position precision
                                            ComputeFloat dx_ov_dt,
                                            ComputeFloat dy_ov_dt,
                                            int          i_domain_begin,
                                            int          j_domain_begin,
                                            int          nprimy,
                                            int          not_spectral_,
                                            unsigned int oversize_,
                                            bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize;
            const unsigned int bin_count      = gridDim.x * gridDim.y;

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y + y_cluster_coordinate;
            const unsigned int thread_index_offset           = threadIdx.x;

            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 2 );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 2 );

            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 2, 2 );
            const int PaddedGCWidth = GPUClusterWithGCWidth + static_cast<int>( kShmemPad ); // To reduce Bank Conflict

            static constexpr unsigned int kLogicalGCWidth   = Params::getGPUClusterWithGhostCellWidth( 2, 2 );
            static constexpr unsigned int kPaddedGCWidth    = kLogicalGCWidth + kShmemPad;
            static constexpr unsigned int kPaddedFieldSize  = kPaddedGCWidth * kPaddedGCWidth;

            // NOTE: I tried having only one cache and reusing it. Doing that
            // requires you to iterate multiple time over the particle which is
            // possible but cost more bandwidth. The speedup was ~x0.92.
            __shared__ ReductionFloat Jx_scratch_space[kNBufferCopies][kPaddedFieldSize]; // Copie to reduce contention
            __shared__ ReductionFloat Jy_scratch_space[kNBufferCopies][kPaddedFieldSize];
            __shared__ ReductionFloat Jz_scratch_space[kNBufferCopies][kPaddedFieldSize];
            __shared__ ReductionFloat rho_scratch_space[kNBufferCopies][kPaddedFieldSize];
            //extern __shared__ ReductionFloat Jrho_scratch_space[];
            //ReductionFloat *Jx_scratch_space = Jrho_scratch_space;
            //ReductionFloat *Jy_scratch_space = (ReductionFloat*)&Jx_scratch_space[kFieldScratchSpaceSize]; //ReductionFloat *Jy_scratch_space = Jrho_scratch_space[1*kFieldScratchSpaceSize];
            //ReductionFloat *Jz_scratch_space = (ReductionFloat*)&Jy_scratch_space[kFieldScratchSpaceSize]; //ReductionFloat *Jz_scratch_space = Jrho_scratch_space[2*kFieldScratchSpaceSize];
            //ReductionFloat *rho_scratch_space = (ReductionFloat*)&Jz_scratch_space[kFieldScratchSpaceSize]; //ReductionFloat *rho_scratch_space = Jrho_scratch_space[3*kFieldScratchSpaceSize];

            const unsigned int buffer_copy_id = threadIdx.x % kNBufferCopies;

            for( unsigned int copy = 0; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
                    Jx_scratch_space[copy][field_index]  = static_cast<ReductionFloat>( 0.0 );
                    Jy_scratch_space[copy][field_index]  = static_cast<ReductionFloat>( 0.0 );
                    Jz_scratch_space[copy][field_index]  = static_cast<ReductionFloat>( 0.0 );
                    rho_scratch_space[copy][field_index] = static_cast<ReductionFloat>( 0.0 );
                }
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
                const ComputeFloat invgf                  = static_cast<ComputeFloat>( device_invgf_[particle_index] );
                const int *const __restrict__ iold        = &device_iold_[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];

                // ── S0 coefficients ──────────────────────────────────
                {
                    const ComputeFloat delta  = static_cast<ComputeFloat>( deltaold[0 * particle_count] );  // explicit cast
                    const ComputeFloat delta2 = delta * delta;
                    Sx0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx0[4] = static_cast<ComputeFloat>( 0.0 );
                }
                {
                    const ComputeFloat delta  = static_cast<ComputeFloat>( deltaold[1 * particle_count] );  // explicit cast
                    const ComputeFloat delta2 = delta * delta;
                    Sy0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy0[4] = static_cast<ComputeFloat>( 0.0 );
                }

                // ── S1 coefficients ──────────────────────────────────
                // [v4] Position calculation in ComputePositionFloat
                {
                    const ComputePositionFloat xpn   = static_cast<ComputePositionFloat>( device_particle_position_x[particle_index] ) * dx_inv;  //
                    const int          ip            = std::round( xpn );                                                                          //
                    const int          ipo           = iold[0 * particle_count];
                    const int          ip_m_ipo      = ip - ipo - i_domain_begin;
                    const ComputeFloat delta         = static_cast<ComputeFloat>( xpn - static_cast<ComputePositionFloat>( ip ) );                 // [v4]
                    const ComputeFloat delta2        = delta * delta;

                    Sx1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[1] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sx1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sx1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sx1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sx1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }
                {
                    const ComputePositionFloat ypn   = static_cast<ComputePositionFloat>( device_particle_position_y[particle_index] ) * dy_inv;  //
                    const int          jp            = std::round( ypn );                                                                          //
                    const int          jpo           = iold[1 * particle_count];
                    const int          jp_m_jpo      = jp - jpo - j_domain_begin;
                    const ComputeFloat delta         = static_cast<ComputeFloat>( ypn - static_cast<ComputePositionFloat>( jp ) );                 //
                    const ComputeFloat delta2        = delta * delta;

                    Sy1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[1] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sy1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sy1[jp_m_jpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sy1[jp_m_jpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sy1[jp_m_jpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }

                // ── Deposition weights ───────────────────────────────
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * static_cast<ComputeFloat>( 1.0 / 3.0 ) * static_cast<ComputeFloat>( device_particle_momentum_z[particle_index] ) * invgf;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */-
                                (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */-
                                (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;

                // ── Jx deposition ────────────────────────────────────
                ComputeFloat tmpJx[5]{};
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( static_cast<ComputeFloat>( 0.5 ) * ( Sy1[0] - Sy0[0] ) );
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( tmpJx[0] ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + static_cast<ComputeFloat>( 0.5 ) * ( Sy1[j] - Sy0[j] ) );
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmpJx[j] ) );
                    }
                }

                // ── Jy deposition ────────────────────────────────────
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int    iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; j++ ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int    iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    ComputeFloat tmp{};
                    for( unsigned int j = 1; j < 5; ++j ) {
                        tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + static_cast<ComputeFloat>( 0.5 ) * ( Sx1[i] - Sx0[i] ) );
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( tmp ) );
                    }
                }

                // ── Jz deposition ────────────────────────────────────
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] ) +
                                                                                                                      Sy1[j] * ( Sx1[i] ) ) ) );
                    }
                }
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( crz_p * ( Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( crz_p * ( Sy0[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx1[i] + Sx0[i] ) +
                                                                                                                      Sy1[j] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i] + Sx1[i] ) ) ) );
                    }
                }

                // ── Rho deposition ───────────────────────────────────
                for( unsigned int i = 0; i < 1; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &rho_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( charge_weight * ( Sx1[0] * Sy1[i] ) ) );
                    for( unsigned int j = 1; j < 5; j++ ) {
                        atomic::LDS::AddNoReturn( &rho_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( charge_weight * ( Sx1[0] * Sy1[j] ) ) );
                    }
                }
                for( unsigned int i = 1; i < 5; ++i ) {
                    const int iloc = ( i + ipo ) * PaddedGCWidth + jpo;
                    atomic::LDS::AddNoReturn( &rho_scratch_space[buffer_copy_id][iloc], static_cast<ReductionFloat>( charge_weight * ( Sx1[i] * Sy1[0] ) ) );
                    for( unsigned int j = 1; j < 5; ++j ) {
                        atomic::LDS::AddNoReturn( &rho_scratch_space[buffer_copy_id][iloc + j], static_cast<ReductionFloat>( charge_weight * ( Sx1[i] * Sy1[j] ) ) );
                    }
                }
            }

            __syncthreads();

            // Reduction
            for( unsigned int copy = 1; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
                    Jx_scratch_space[0][field_index]  += Jx_scratch_space[copy][field_index];
                    Jy_scratch_space[0][field_index]  += Jy_scratch_space[copy][field_index];
                    Jz_scratch_space[0][field_index]  += Jz_scratch_space[copy][field_index];
                    rho_scratch_space[0][field_index] += rho_scratch_space[copy][field_index];
                }
            }

            if( kNBufferCopies > 1 ) {
                __syncthreads();
            }

            // Flush to global memory
            for( unsigned int field_index = thread_index_offset;
                 field_index < kPaddedFieldSize;
                 field_index += workgroup_size ) {

                const unsigned int local_x = field_index / PaddedGCWidth;
                const unsigned int local_y = field_index % PaddedGCWidth;

                if( kShmemPad > 0 ) {
                    if( local_x >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_y >= static_cast<unsigned int>(GPUClusterWithGCWidth) ) {
                        continue;
                    }
                }

                const unsigned int global_x = global_x_scratch_space_coordinate_offset + local_x + (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;
                const unsigned int global_y = global_y_scratch_space_coordinate_offset + local_y + (oversize_-2) /* Offset to get in-cluster coordinate with an oversize > 2*/ ;

                const unsigned int global_memory_index = global_x * nprimy + global_y;

                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + not_spectral_ * global_x], static_cast<double>( Jy_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index], static_cast<double>( Jz_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_rho[global_memory_index], static_cast<double>( rho_scratch_space[0][field_index] ) );
            }
        }
    } // namespace kernel


    // ====================================================================
    // Launch wrapper: Current deposition (Jx, Jy, Jz)
    // ====================================================================

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
                             int    not_spectral_,
                             unsigned int oversize_,
                             bool   cell_sorting )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 2 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), 1 };

        // ================================================================
        // ▼▼▼ TUNING PARAMETERS — ADJUST THESE ▼▼▼
        // ================================================================

        static constexpr std::size_t kWorkgroupSize   = 256;
        static constexpr std::size_t kNBufferCopies   = 1; // Very important to reduce contention, 1: Legacy. 8: is good
        static constexpr std::size_t kShmemPad        = 0; // 0: Legacy. 2: NCU give a local minimum
        static constexpr std::size_t kMinBlocksPerSM  = 1; // 1: Legacy. 6: NCU give a local minimum

        //static constexpr std::size_t kWorkgroupSize   = 256;
        //static constexpr std::size_t kNBufferCopies   = 8;
        //static constexpr std::size_t kShmemPad        = 2;
        //static constexpr std::size_t kMinBlocksPerSM  = 6;

        // ================================================================
        // ▲▲▲ END TUNING PARAMETERS ▲▲▲
        // ================================================================

        const ::dim3 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // Type aliases — change ComputeFloat to float for mixed precision
        using ComputeFloat         = double;  // ← set to float for mixed precision mode
        using ComputePositionFloat = double;  // ← always double (position safety)
        using ReductionFloat       = double;  // ← set to float for mixed precision mode

        auto KernelFunction = kernel::DepositCurrentDensity_2D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
#if defined ( __HIP__ ) 
        hipLaunchKernelGGL( KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, 0,
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
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),    // kept as double
                            static_cast<ComputePositionFloat>( dy_inv ),    // kept as double
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_,
                            oversize_,
                            cell_sorting );
        checkHIPErrors( ::hipDeviceSynchronize() );
#elif defined (  __NVCC__ )
        KernelFunction <<<
                            kGridDimension,
                            kBlockDimension,
                            0, 0
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
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),    // kept as double
                            static_cast<ComputePositionFloat>( dy_inv ),    // kept as double
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif
    }


    // ====================================================================
    // Launch wrapper: Current + Density deposition (Jx, Jy, Jz, rho)
    // ====================================================================

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
                                       int    not_spectral_,
                                       unsigned int oversize_,
                                       bool cell_sorting )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 2 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        const ::dim3 kGridDimension  { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), 1 };

        // ================================================================
        // ▼▼▼ TUNING PARAMETERS — ADJUST THESE ▼▼▼
        // ================================================================

        static constexpr std::size_t kWorkgroupSize   = 256;
        static constexpr std::size_t kNBufferCopies   = 1; // Very important to reduce contention, 1: Legacy. 8: is good
        static constexpr std::size_t kShmemPad        = 0; // 0: Legacy. 2: NCU give a local minimum
        static constexpr std::size_t kMinBlocksPerSM  = 1; // 1: Legacy. 6: NCU give a local minimum

        //static constexpr std::size_t kWorkgroupSize   = 256;
        //static constexpr std::size_t kNBufferCopies   = 8;
        //static constexpr std::size_t kShmemPad        = 2;
        //static constexpr std::size_t kMinBlocksPerSM  = 6;

        // ================================================================
        // ▲▲▲ END TUNING PARAMETERS ▲▲▲
        // ================================================================

        const ::dim3 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // Type aliases — change ComputeFloat to float for mixed precision
        using ComputeFloat         = double;  // ← set to float for mixed precision mode
        using ComputePositionFloat = double;  // ← always double (position safety)
        using ReductionFloat       = double;  // ← set to float for mixed precision mode

        auto KernelFunction = kernel::DepositCurrentAndDensity_2D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
#if defined ( __HIP__ ) 
        hipLaunchKernelGGL( KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, 0,
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
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),    // kept as double
                            static_cast<ComputePositionFloat>( dy_inv ),    // kept as double
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_,
                            oversize_,
                            cell_sorting );
        checkHIPErrors( ::hipDeviceSynchronize() );
#elif defined (  __NVCC__ )
        KernelFunction <<<
                            kGridDimension,
                            kBlockDimension,
                            0, 0
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
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),    // kept as double
                            static_cast<ComputePositionFloat>( dy_inv ),    // kept as double
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            i_domain_begin, j_domain_begin,
                            nprimy,
                            not_spectral_,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif 
    }

} // namespace cudahip2d

#endif
