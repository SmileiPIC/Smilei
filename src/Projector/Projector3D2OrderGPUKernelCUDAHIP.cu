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
//!   kNBufferCopies   — Replicated shared memory buffers (1 = off) to reduce contention
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

#if defined( SMILEI_ACCELERATOR_GPU )

#if defined( __HIP__ )
    #include <hip/hip_runtime.h>
#elif defined( __NVCC__ )
    #include <cuda_runtime.h>
    #include <cuda.h>
#endif

#include "Params.h"
#include "gpu.h"
#include "stdio.h"
#include <iostream>
namespace cudahip {
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

        // ── Helper: S0 shape function ────────────────────────────
        // deltaold is already a delta in [-0.5, 0.5], safe to use in ComputeFloat.
        // Takes double deltaold_value from memory, casts to ComputeFloat.
        template <typename ComputeFloat>
        __device__ void inline __attribute__((always_inline)) init_S0(const double deltaold_value, ComputeFloat *__restrict__ S0)
        {
            const ComputeFloat delta  = static_cast<ComputeFloat>( deltaold_value );  // explicit cast from double
            const ComputeFloat delta2 = delta * delta;
            S0[0] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
            S0[1] = static_cast<ComputeFloat>( 0.75 ) - delta2;
            S0[2] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
            S0[3] = static_cast<ComputeFloat>( 0.0 ) ;
        }

        // ── Helper: S1 shape function ────────────────────────────
        // CRITICAL: xpn is computed in ComputePositionFloat (double) to
        // preserve precision during the catastrophic cancellation xpn - ip.
        // The resulting delta is then narrowed to ComputeFloat for shapes.
        template <typename ComputeFloat, typename ComputePositionFloat>
        __device__ void inline __attribute__((always_inline)) init_S1(const ComputePositionFloat xpn, const int ipo, const int i_domain_begin,
                                                                      ComputeFloat *__restrict__ S1)
        {
            const int          ip       = std::round( xpn );                                                       // round in high precision
            const int          ip_m_ipo = ip - ipo - i_domain_begin;
            const ComputeFloat delta    = static_cast<ComputeFloat>( xpn - static_cast<ComputePositionFloat>( ip ) ); // subtract in high precision, then narrow
            const ComputeFloat delta2   = delta * delta;

            S1[0] = static_cast<ComputeFloat>( 0.0 );
            S1[1] = static_cast<ComputeFloat>( 0.0 );
            S1[3] = static_cast<ComputeFloat>( 0.0 );
            S1[4] = static_cast<ComputeFloat>( 0.0 );

            S1[ip_m_ipo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
            S1[ip_m_ipo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
            S1[ip_m_ipo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
        }


        // ====================================================================
        // Current deposition kernel (Jx, Jy, Jz)
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
        DepositCurrentDensity_3D_Order2( double *__restrict__ device_Jx,
                                         double *__restrict__ device_Jy,
                                         double *__restrict__ device_Jz,
                                         int Jx_size,
                                         int Jy_size,
                                         int Jz_size,
                                         const double *__restrict__ device_particle_position_x,
                                         const double *__restrict__ device_particle_position_y,
                                         const double *__restrict__ device_particle_position_z,
                                         const short *__restrict__ device_particle_charge,
                                         const double *__restrict__ device_particle_weight,
                                         const int *__restrict__ device_bin_index,
                                         const double *__restrict__ device_invgf_,
                                         int *__restrict__ device_iold,
                                         const double *__restrict__ device_deltaold_,
                                         ComputeFloat inv_cell_volume,
                                         ComputePositionFloat dx_inv,  // position precision
                                         ComputePositionFloat dy_inv,  // position precision
                                         ComputePositionFloat dz_inv,  // position precision
                                         ComputeFloat dx_ov_dt,
                                         ComputeFloat dy_ov_dt,
                                         ComputeFloat dz_ov_dt,
                                         int          i_domain_begin,
                                         int          j_domain_begin,
                                         int          k_domain_begin,
                                         int          nprimy,
                                         int          nprimz,
                                         int          not_spectral,
                                         unsigned int oversize_,
                                         bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize;
            const unsigned int bin_count      = gridDim.x * gridDim.y * gridDim.z;

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int z_cluster_coordinate          = blockIdx.z;

            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y * gridDim.z + y_cluster_coordinate * gridDim.z + z_cluster_coordinate;
            const unsigned int thread_index_offset           = threadIdx.x;

            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 3 );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 3 );
            const unsigned int global_z_scratch_space_coordinate_offset = z_cluster_coordinate * Params::getGPUClusterWidth( 3 );

            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 3, 2 );
            const int PaddedGCWidth = GPUClusterWithGCWidth + static_cast<int>( kShmemPad ); // To reduce Bank Conflict

            ComputeFloat one_third = 1. / 3.;

            // NOTE: We gain from the particles not being sorted inside a
            // cluster because it reduces the bank conflicts one gets when
            // multiple threads access the same part of the shared memory. Such
            // "conflicted" accesses are serialized !
            // NOTE: We use a bit to much LDS. For Jx, the first row could be
            // discarded, for Jy we could remove the first column.

            static constexpr unsigned int kLogicalGCWidth    = Params::getGPUClusterWithGhostCellWidth( 3, 2 );
            static constexpr unsigned int kPaddedGCWidth     = kLogicalGCWidth + kShmemPad;
            static constexpr unsigned int kPaddedFieldSize   = kPaddedGCWidth * kPaddedGCWidth * kPaddedGCWidth;
            static constexpr unsigned int kLogicalFieldSize  = kLogicalGCWidth * kLogicalGCWidth * kLogicalGCWidth;

            __shared__ ReductionFloat Jx_scratch_space[kNBufferCopies][kPaddedFieldSize]; // Copie to reduce contention
            __shared__ ReductionFloat Jy_scratch_space[kNBufferCopies][kPaddedFieldSize];
            __shared__ ReductionFloat Jz_scratch_space[kNBufferCopies][kPaddedFieldSize];

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

            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 :
                                                                                     device_bin_index[workgroup_dedicated_bin_index - 1];
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
                const int *const __restrict__ iold        = &device_iold[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[4];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[4];
                ComputeFloat Sy1[5];
                ComputeFloat Sz0[4];
                ComputeFloat Sz1[5];

                // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

                // S0: deltaold read as double, cast to ComputeFloat inside init_S0
                init_S0( deltaold[0 * particle_count], Sx0 );
                init_S0( deltaold[1 * particle_count], Sy0 );
                init_S0( deltaold[2 * particle_count], Sz0 );

                // S1: position × d_inv computed in ComputePositionFloat, delta narrowed to ComputeFloat
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_x[particle_index] ) * dx_inv,
                         iold[0 * particle_count], i_domain_begin, Sx1 );
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_y[particle_index] ) * dy_inv,
                         iold[1 * particle_count], j_domain_begin, Sy1 );
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_z[particle_index] ) * dz_inv,
                         iold[2 * particle_count], k_domain_begin, Sz1 );

                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * dz_ov_dt;

             
                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);
                const int kpo = iold[2 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_z_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);
 
                // ========================================================
                // Jx deposition
                // ========================================================

                // j=0, k=0
                {
                    ComputeFloat tmp = crx_p * Sy1[0] * one_third * Sz1[0];
                    ComputeFloat tmp_reduction{};
                    const int jk_loc = ( ipo * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo;
                    tmp_reduction -= Sx1[0] * tmp;
                    const int loc = PaddedGCWidth*PaddedGCWidth + jk_loc;
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );

                    for( unsigned int i = 2; i < 5; ++i ) {
                        tmp_reduction -= ( Sx1[i-1] - Sx0[i-2] ) * tmp;
                        const int loc = i*PaddedGCWidth*PaddedGCWidth + jk_loc;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // j=0, k=1..3
                for( unsigned int k = 1; k < 4; ++k ) {
                    ComputeFloat tmp = crx_p * Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sz0[k-1]
                                                    + one_third * ( Sz1[k] - Sz0[k-1] ) );
                    ComputeFloat tmp_reduction{};
                    const int jk_loc = ( ipo * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo + k;
                    tmp_reduction -= Sx1[0] * tmp;
                    const int loc = PaddedGCWidth*PaddedGCWidth + jk_loc;
                    atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );

                    for( unsigned int i = 2; i < 5; ++i ) {
                        tmp_reduction -= ( Sx1[i-1] - Sx0[i-2] ) * tmp;
                        const int loc = i*PaddedGCWidth*PaddedGCWidth + jk_loc;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // j=1..3
                for( unsigned int j = 1; j < 4; ++j ) {
                    // k=0
                    {
                        ComputeFloat tmp = crx_p * ( static_cast<ComputeFloat>( 0.5 ) * Sz1[0]*Sy0[j-1] 
                                                     + one_third * ( Sy1[j] - Sy0[j-1] ) * Sz1[0] );
                        ComputeFloat tmp_reduction{};
                        const int jk_loc = ( ipo * PaddedGCWidth + jpo + j ) * PaddedGCWidth + kpo;
                        tmp_reduction -= Sx1[0] * tmp;
                        const int loc = PaddedGCWidth*PaddedGCWidth + jk_loc;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int i = 2; i < 5; ++i ) {
                            tmp_reduction -= ( Sx1[i-1] - Sx0[i-2] ) * tmp;
                            const int loc = i*PaddedGCWidth*PaddedGCWidth + jk_loc;
                            atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                    // k=1..3
                    for( unsigned int k = 1; k < 4; ++k ) {
                        ComputeFloat tmp = crx_p * (   Sy0[j-1]*Sz0[k-1]
                                                     + static_cast<ComputeFloat>( 0.5 ) * ( ( Sy1[j] - Sy0[j-1] )*Sz0[k-1] + ( Sz1[k] - Sz0[k-1] )*Sy0[j-1] )
                                                     + one_third * ( Sy1[j] - Sy0[j-1] ) * ( Sz1[k] - Sz0[k-1] ) );
                        ComputeFloat tmp_reduction{};
                        const int jk_loc = ( ipo * PaddedGCWidth + jpo + j ) * PaddedGCWidth + kpo + k;
                        tmp_reduction -= Sx1[0] * tmp;
                        const int loc = PaddedGCWidth*PaddedGCWidth + jk_loc;
                        atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int i = 2; i < 5; ++i ) {
                            tmp_reduction -= ( Sx1[i-1] - Sx0[i-2] ) * tmp;
                            const int loc = i*PaddedGCWidth*PaddedGCWidth + jk_loc;
                            atomic::LDS::AddNoReturn( &Jx_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // ========================================================
                // Jy deposition
                // ========================================================

                // i=0, k=0
                {
                    ComputeFloat tmp = cry_p * Sx1[0] * one_third * Sz1[0];
                    ComputeFloat tmp_reduction{};
                    const int ik_loc = (( 0 + ipo ) * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo;
                    tmp_reduction -= Sy1[0] * tmp;
                    const int loc = PaddedGCWidth + ik_loc;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int j = 2; j < 5; ++j ) {
                        tmp_reduction -= ( Sy1[j-1] - Sy0[j-2] ) * tmp;
                        const int loc = j*PaddedGCWidth + ik_loc;
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // i=0, k=1..3
                for( unsigned int k = 1; k < 4; ++k ) {
                    ComputeFloat tmp = cry_p * Sx1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sz0[k-1]
                                                    + one_third * ( Sz1[k] - Sz0[k-1] ) );
                    ComputeFloat tmp_reduction{};
                    const int ik_loc = (( 0 + ipo ) * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo + k;
                    tmp_reduction -= Sy1[0] * tmp;
                    const int loc = PaddedGCWidth + ik_loc;
                    atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int j = 2; j < 5; ++j ) {
                        tmp_reduction -= ( Sy1[j-1] - Sy0[j-2] ) * tmp;
                        const int loc = j*PaddedGCWidth + ik_loc;
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // i=1..3
                for( unsigned int i = 1; i < 4; ++i ) {
                    // k=0
                    {
                        ComputeFloat tmp = cry_p * ( static_cast<ComputeFloat>( 0.5 ) * Sz1[0] * Sx0[i-1] 
                                                     + one_third * ( Sx1[i] - Sx0[i-1] ) * Sz1[0] );
                        ComputeFloat tmp_reduction{};
                        const int ik_loc = (( i + ipo ) * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo + 0;
                        tmp_reduction -= Sy1[0] * tmp;
                        const int loc = PaddedGCWidth + ik_loc;
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int j = 2; j < 5; ++j ) {
                            tmp_reduction -= ( Sy1[j-1] - Sy0[j-2] ) * tmp;
                            const int loc = j*PaddedGCWidth + ik_loc;
                            atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                    // k=1..3
                    for( unsigned int k = 1; k < 4; ++k ) {
                        ComputeFloat tmp = cry_p * (   Sx0[i-1]*Sz0[k-1]
                                                     + static_cast<ComputeFloat>( 0.5 ) * ( ( Sx1[i] - Sx0[i-1] )*Sz0[k-1] + ( Sz1[k] - Sz0[k-1] )*Sx0[i-1] )
                                                     + one_third * ( Sx1[i] - Sx0[i-1] ) * ( Sz1[k] - Sz0[k-1] ) );
                        ComputeFloat tmp_reduction{};
                        const int ik_loc = (( i + ipo ) * PaddedGCWidth + jpo ) * PaddedGCWidth + kpo + k;
                        tmp_reduction -= Sy1[0] * tmp;
                        const int loc = PaddedGCWidth + ik_loc;
                        atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int j = 2; j < 5; ++j ) {
                            tmp_reduction -= ( Sy1[j-1] - Sy0[j-2] ) * tmp;
                            const int loc = j*PaddedGCWidth + ik_loc;
                            atomic::LDS::AddNoReturn( &Jy_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // ========================================================
                // Jz deposition
                // ========================================================

                // i=0, j=0
                {
                    ComputeFloat tmp = crz_p * one_third * Sx1[0] * Sy1[0];
                    ComputeFloat tmp_reduction{};
                    const int ij_loc = (( 0 + ipo ) * PaddedGCWidth + (jpo + 0 )) * PaddedGCWidth + kpo;
                    tmp_reduction -= Sz1[0] * tmp;
                    const int loc = 1 + ij_loc;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int k = 2; k < 5; ++k ) {
                        tmp_reduction -= ( Sz1[k-1] - Sz0[k-2] ) * tmp;
                        const int loc = k + ij_loc;
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // i=0, j=1..3
                for( unsigned int j = 1; j < 4; ++j ) {
                    ComputeFloat tmp = crz_p * Sx1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sy0[j-1]
                                                    + one_third * ( Sy1[j] - Sy0[j-1] ) );
                    ComputeFloat tmp_reduction{};
                    const int ij_loc = (( 0 + ipo ) * PaddedGCWidth + (jpo + j)) * PaddedGCWidth + kpo;
                    tmp_reduction -= Sz1[0] * tmp;
                    const int loc = 1 + ij_loc;
                    atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    for( unsigned int k = 2; k < 5; ++k ) {
                        tmp_reduction -= ( Sz1[k-1] - Sz0[k-2] ) * tmp;
                        const int loc = k + ij_loc;
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                    }
                }
                // i=1..3
                for( unsigned int i = 1; i < 4; ++i ) {
                    // j=0
                    {
                        ComputeFloat tmp = crz_p * Sy1[0] * ( static_cast<ComputeFloat>( 0.5 ) * Sx0[i-1]
                                                     + one_third * ( Sx1[i] - Sx0[i-1] ) );
                        ComputeFloat tmp_reduction{};
                        const int ij_loc = (( i + ipo ) * PaddedGCWidth + (jpo + 0 )) * PaddedGCWidth + kpo;
                        tmp_reduction -= Sz1[0] * tmp;
                        const int loc = 1 + ij_loc;
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int k = 2; k < 5; ++k ) {
                            tmp_reduction -= ( Sz1[k-1] - Sz0[k-2] ) * tmp;
                            const int loc = k + ij_loc;
                            atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                    // j=1..3
                    for( unsigned int j = 1; j < 4; ++j ) {
                        ComputeFloat tmp = crz_p * (   Sx0[i-1]*Sy0[j-1]
                                                     + static_cast<ComputeFloat>( 0.5 ) * ( ( Sx1[i] - Sx0[i-1] )*Sy0[j-1] + ( Sy1[j] - Sy0[j-1] )*Sx0[i-1] )
                                                     + one_third * ( Sx1[i] - Sx0[i-1] ) * ( Sy1[j] - Sy0[j-1] ) );
                        ComputeFloat tmp_reduction{};
                        const int ij_loc = (( i + ipo ) * PaddedGCWidth + (jpo + j)) * PaddedGCWidth + kpo;
                        tmp_reduction -= Sz1[0] * tmp;
                        const int loc = 1 + ij_loc;
                        atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        for( unsigned int k = 2; k < 5; ++k ) {
                            tmp_reduction -= ( Sz1[k-1] - Sz0[k-2] ) * tmp;
                            const int loc = k + ij_loc;
                            atomic::LDS::AddNoReturn( &Jz_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }
            } // end particle loop

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

                const unsigned int local_x = field_index / ( PaddedGCWidth * PaddedGCWidth );
                const unsigned int local_y = ( field_index % ( PaddedGCWidth * PaddedGCWidth ) ) / PaddedGCWidth;
                const unsigned int local_z = field_index % PaddedGCWidth;

                if( kShmemPad > 0 ) {
                    if( local_x >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_y >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_z >= static_cast<unsigned int>(GPUClusterWithGCWidth) ) {
                        continue;
                    }
                }

                const unsigned int global_x = global_x_scratch_space_coordinate_offset + local_x +(oversize_-2);
                const unsigned int global_y = global_y_scratch_space_coordinate_offset + local_y +(oversize_-2);
                const unsigned int global_z = global_z_scratch_space_coordinate_offset + local_z +(oversize_-2);

                const unsigned int global_memory_index = ( global_x * nprimy + global_y ) * nprimz + global_z;

                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index],
                                          static_cast<double>( Jx_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + not_spectral * global_x * nprimz],
                                          static_cast<double>( Jy_scratch_space[0][field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index + not_spectral * (global_x * nprimy + global_y)],
                                          static_cast<double>( Jz_scratch_space[0][field_index] ) );
            }
        } // end DepositCurrentDensity


        // ====================================================================
        // Density deposition kernel (rho)
        // ====================================================================

        template <typename ComputeFloat,
                  typename ComputePositionFloat,       // [v4] NEW
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize,
                  std::size_t kNBufferCopies,
                  std::size_t kShmemPad,
                  std::size_t kMinBlocksPerSM>
        __global__ void
        __launch_bounds__(kWorkgroupSize, kMinBlocksPerSM)
        DepositDensity_3D_Order2(
                                            double *__restrict__ device_rho,
                                            int rho_size,
                                            const double *__restrict__ device_particle_position_x,
                                            const double *__restrict__ device_particle_position_y,
                                            const double *__restrict__ device_particle_position_z,
                                            const short *__restrict__ device_particle_charge,
                                            const double *__restrict__ device_particle_weight,
                                            const int *__restrict__ device_bin_index,
                                            const double *__restrict__ device_invgf_,
                                            int *__restrict__ device_iold,
                                            const double *__restrict__ device_deltaold_,
                                            ComputeFloat inv_cell_volume,
                                            ComputePositionFloat dx_inv,  // [v4] position precision
                                            ComputePositionFloat dy_inv,  // [v4] position precision
                                            ComputePositionFloat dz_inv,  // [v4] position precision
                                            ComputeFloat dx_ov_dt,
                                            ComputeFloat dy_ov_dt,
                                            ComputeFloat dz_ov_dt,
                                            int          i_domain_begin,
                                            int          j_domain_begin,
                                            int          k_domain_begin,
                                            int          nprimy,
                                            int          nprimz,
                                            int          not_spectral,
                                            unsigned int oversize_,
                                            bool         cell_sorting )
        {
            const unsigned int workgroup_size = kWorkgroupSize;
            const unsigned int bin_count      = gridDim.x * gridDim.y * gridDim.z;

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int z_cluster_coordinate          = blockIdx.z;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y * gridDim.z + y_cluster_coordinate * gridDim.z + z_cluster_coordinate;
            const unsigned int thread_index_offset           = threadIdx.x;

            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 3 );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 3 );
            const unsigned int global_z_scratch_space_coordinate_offset = z_cluster_coordinate * Params::getGPUClusterWidth( 3 );

            const int GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 3, 2 );
            const int PaddedGCWidth = GPUClusterWithGCWidth + static_cast<int>( kShmemPad ); // Reduce Bank Conflict
            ComputeFloat one_third = 1. / 3.;

            static constexpr unsigned int kLogicalGCWidth    = Params::getGPUClusterWithGhostCellWidth( 3, 2 );
            static constexpr unsigned int kPaddedGCWidth     = kLogicalGCWidth + kShmemPad;
            static constexpr unsigned int kPaddedFieldSize   = kPaddedGCWidth * kPaddedGCWidth * kPaddedGCWidth;

            __shared__ ReductionFloat rho_scratch_space[kNBufferCopies][kPaddedFieldSize]; // Copie to reduce contention

            const unsigned int buffer_copy_id = threadIdx.x % kNBufferCopies;

            for( unsigned int copy = 0; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
                    rho_scratch_space[copy][field_index] = static_cast<ReductionFloat>( 0.0 );
                }
            }

            __syncthreads();

            const unsigned int particle_count = device_bin_index[bin_count - 1];

            const unsigned int first_particle = workgroup_dedicated_bin_index == 0 ? 0 :
                                                                                     device_bin_index[workgroup_dedicated_bin_index - 1];
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
                const int *const __restrict__ iold        = &device_iold[particle_index];

                ComputeFloat Sx1[5];
                ComputeFloat Sy1[5];
                ComputeFloat Sz1[5];

                // [v4] S1: position × d_inv computed in ComputePositionFloat
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_x[particle_index] ) * dx_inv,
                         iold[0 * particle_count], i_domain_begin, Sx1 );
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_y[particle_index] ) * dy_inv,
                         iold[1 * particle_count], j_domain_begin, Sy1 );
                init_S1( static_cast<ComputePositionFloat>( device_particle_position_z[particle_index] ) * dz_inv,
                         iold[2 * particle_count], k_domain_begin, Sz1 );

                const ComputeFloat charge_weight = inv_cell_volume * 
                                                   static_cast<ComputeFloat>( device_particle_charge[particle_index] ) *
                                                   static_cast<ComputeFloat>( device_particle_weight[particle_index] );

                // const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                // const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                // const ComputeFloat crz_p         = charge_weight * dz_ov_dt;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);
                const int kpo = iold[2 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_z_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */ -
                                (oversize_-2);

                // Rho deposition
                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int j = 0; j < 5; ++j ) {
		                for( unsigned int k = 0; k < 5; ++k ) {
			                ComputeFloat tmp = charge_weight * Sx1[i]*Sy1[j];
			                const int ij_loc = (( i + ipo ) * PaddedGCWidth +
                                (jpo + j)) * PaddedGCWidth + kpo; 
                            const int loc = ij_loc + k;
                            atomic::LDS::AddNoReturn( &rho_scratch_space[buffer_copy_id][loc], static_cast<ReductionFloat>( tmp * Sz1[k] ) );
                        }
                    }
                }
            } // end particle loop

            __syncthreads();

            // Reduction (no-op when kNBufferCopies=1)
            for( unsigned int copy = 1; copy < kNBufferCopies; ++copy ) {
                for( unsigned int field_index = thread_index_offset;
                     field_index < kPaddedFieldSize;
                     field_index += workgroup_size ) {
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

                const unsigned int local_x = field_index / ( PaddedGCWidth * PaddedGCWidth );
                const unsigned int local_y = ( field_index % ( PaddedGCWidth * PaddedGCWidth ) ) / PaddedGCWidth;
                const unsigned int local_z = field_index % PaddedGCWidth;

                if( kShmemPad > 0 ) {
                    if( local_x >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_y >= static_cast<unsigned int>(GPUClusterWithGCWidth) ||
                        local_z >= static_cast<unsigned int>(GPUClusterWithGCWidth) ) {
                        continue;
                    }
                }

                const unsigned int global_x = global_x_scratch_space_coordinate_offset + local_x + (oversize_-2) ;
                const unsigned int global_y = global_y_scratch_space_coordinate_offset + local_y + (oversize_-2) ;
                const unsigned int global_z = global_z_scratch_space_coordinate_offset + local_z + (oversize_-2) ;

                const unsigned int global_memory_index = ( global_x * nprimy + global_y ) * nprimz + global_z;

                atomic::GDS::AddNoReturn( &device_rho[global_memory_index], static_cast<double>( rho_scratch_space[0][field_index] ) );
            }
        }


   } // namespace kernel


    // ====================================================================
    // Launch wrapper: Current deposition (Jx, Jy, Jz)
    // ====================================================================

    void
    currentDepositionKernel3D( double *__restrict__ host_Jx,
                               double *__restrict__ host_Jy,
                               double *__restrict__ host_Jz,
                               int Jx_size,
                               int Jy_size,
                               int Jz_size,
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
                               int    not_spectral,
                               unsigned int oversize_,
                               bool   cell_sorting )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 3 /* 3D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        const ::dim3 kGridDimension { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), static_cast<uint32_t>( z_dimension_bin_count ) };

        // ================================================================
        // ▼▼▼ TUNING PARAMETERS — ADJUST THESE ▼▼▼
        // ================================================================

        // Legacy parameters
        static constexpr std::size_t kWorkgroupSize   = 128;
        static constexpr std::size_t kNBufferCopies   = 1;
        static constexpr std::size_t kShmemPad        = 0;
        static constexpr std::size_t kMinBlocksPerSM  = 1;

        //static constexpr std::size_t kWorkgroupSize   = 128;
        //static constexpr std::size_t kNBufferCopies   = 1;
        //static constexpr std::size_t kShmemPad        = 2;
        //static constexpr std::size_t kMinBlocksPerSM  = 10;

        // ================================================================
        // ▲▲▲ END TUNING PARAMETERS ▲▲▲
        // ================================================================

        const ::dim3 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // Type aliases — change ComputeFloat to float for mixed precision
        using ComputeFloat         = double;  // ← set to float for mixed precision mode
        using ComputePositionFloat = double;  // ← always double (position safety)
        using ReductionFloat       = double;  // ← set to float for mixed precision mode

#if defined ( __HIP__ )
        auto KernelFunction = kernel::DepositCurrentDensity_3D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
        hipLaunchKernelGGL
                        (   KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, 0,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jx ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jy ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_Jz ),
                            Jx_size, Jy_size, Jz_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),
                            static_cast<ComputePositionFloat>( dy_inv ),
                            static_cast<ComputePositionFloat>( dz_inv ),
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            static_cast<ComputeFloat>( dz_ov_dt ),
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral,
                            oversize_,
                            cell_sorting 
                        );
        checkHIPErrors( ::hipDeviceSynchronize() );

#elif defined ( __NVCC__ )
        auto KernelFunction = kernel::DepositCurrentDensity_3D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
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
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),
                            static_cast<ComputePositionFloat>( dy_inv ),
                            static_cast<ComputePositionFloat>( dz_inv ),
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            static_cast<ComputeFloat>( dz_ov_dt ),
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif
    }


    // ====================================================================
    // Launch wrapper: Density deposition (rho)
    // ====================================================================

    void
    densityDepositionKernel3D( 
                                double *__restrict__ host_rho,
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
                                int    not_spectral,
                                unsigned int oversize_,
                                bool   cell_sorting )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 3 /* 3D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        const ::dim3 kGridDimension { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), static_cast<uint32_t>( z_dimension_bin_count ) };

        // ================================================================
        // ▼▼▼ TUNING PARAMETERS — ADJUST THESE ▼▼▼
        // ================================================================

        // Legacy parameters
        static constexpr std::size_t kWorkgroupSize   = 128;
        static constexpr std::size_t kNBufferCopies   = 1;
        static constexpr std::size_t kShmemPad        = 0;
        static constexpr std::size_t kMinBlocksPerSM  = 1;

        //static constexpr std::size_t kWorkgroupSize   = 128;
        //static constexpr std::size_t kNBufferCopies   = 1;
        //static constexpr std::size_t kShmemPad        = 2;
        //static constexpr std::size_t kMinBlocksPerSM  = 10;

        // ================================================================
        // ▲▲▲ END TUNING PARAMETERS ▲▲▲
        // ================================================================

        const ::dim3 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // Type aliases
        using ComputeFloat         = double;   // ← set to float for mixed precision mode
        using ComputePositionFloat = double;   // ← always double (position safety)
        using ReductionFloat       = double;   // ← set to float for mixed precision mode

#if defined ( __HIP__ )
        auto KernelFunction = kernel::DepositDensity_3D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
        
        hipLaunchKernelGGL( KernelFunction,
                            kGridDimension,
                            kBlockDimension,
                            0, 0,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_rho ),
                            rho_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),
                            static_cast<ComputePositionFloat>( dy_inv ),
                            static_cast<ComputePositionFloat>( dz_inv ),
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            static_cast<ComputeFloat>( dz_ov_dt ),
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral,
                            oversize_,
                            cell_sorting );

        checkHIPErrors( ::hipDeviceSynchronize() );
#elif defined (  __NVCC__ )
        auto KernelFunction = kernel::DepositDensity_3D_Order2<ComputeFloat, ComputePositionFloat, ReductionFloat, kWorkgroupSize, kNBufferCopies, kShmemPad, kMinBlocksPerSM>;
        KernelFunction <<<
                            kGridDimension,
                            kBlockDimension,
                            0, 0
                       >>>
                       (
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_rho ),
                            rho_size,
                            device_particle_position_x,
                            device_particle_position_y,
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            static_cast<ComputeFloat>( inv_cell_volume ),
                            static_cast<ComputePositionFloat>( dx_inv ),
                            static_cast<ComputePositionFloat>( dy_inv ),
                            static_cast<ComputePositionFloat>( dz_inv ),
                            static_cast<ComputeFloat>( dx_ov_dt ),
                            static_cast<ComputeFloat>( dy_ov_dt ),
                            static_cast<ComputeFloat>( dz_ov_dt ),
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral,
                            oversize_,
                            cell_sorting
                       );
        checkHIPErrors( ::cudaDeviceSynchronize() );
#endif
    }

} // namespace cudahip

#endif
