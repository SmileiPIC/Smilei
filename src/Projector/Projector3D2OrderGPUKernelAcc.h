//! Optimized Acc projection (from Julien Derouillat) 

#if defined( SMILEI_ACCELERATOR_GPU )

#include <cmath>
#include "Tools.h"
#include "gpu.h"

namespace acc {

    static inline void
    currentDepositionKernel3D( double *__restrict__ Jx,
                               double *__restrict__ Jy,
                               double *__restrict__ Jz,
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
                               const double *__restrict__ invgf,
                               int    *__restrict__ iold,
                               const double *__restrict__ deltaold,
                               const unsigned int    number_of_particles,
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

    const unsigned int bin_count      = 1;
    const int          nparts = host_bin_index[bin_count - 1];

    // TODO(Etienne M): Implement a cuda/hip kernel and enable particle 3D sorting/binning

    const double *const __restrict__ position_x = device_particle_position_x;
    const double *const __restrict__ position_y = device_particle_position_y;
    const double *const __restrict__ position_z = device_particle_position_z;
    const short  *const __restrict__ charge     = device_particle_charge;
    const double *const __restrict__ weight     = device_particle_weight;

    const double one_third = 1./3.;

    const istart = 0; 
    const iend   = nparts; 

    const int packsize = nparts;
    const int npack    = ( ( iend - istart ) + ( packsize - 1 ) ) / packsize; // divide + ceil npack.

    static constexpr bool kAutoDeviceFree = true;
    const std::size_t     kTmpArraySize   = 5 * packsize;

    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sx0{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sy0{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sz0{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sx1{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sy1{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sz1{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_DSx{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_DSy{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_DSz{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_sumX{ kTmpArraySize };

    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sx0 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sy0 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sz0 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sx1 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sy1 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sz1 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_DSx );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_DSy );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_DSz );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_sumX );

    double *const __restrict__ Sx0  = host_device_Sx0.data();
    double *const __restrict__ Sy0  = host_device_Sy0.data();
    double *const __restrict__ Sz0  = host_device_Sz0.data();
    double *const __restrict__ Sx1  = host_device_Sx1.data();
    double *const __restrict__ Sy1  = host_device_Sy1.data();
    double *const __restrict__ Sz1  = host_device_Sz1.data();
    double *const __restrict__ DSx  = host_device_DSx.data();
    double *const __restrict__ DSy  = host_device_DSy.data();
    double *const __restrict__ DSz  = host_device_DSz.data();
    double *const __restrict__ sumX = host_device_sumX.data();

    for (int ipack=0 ; ipack<npack ; ipack++) {
        const int istart_pack       = istart + ipack * packsize;
        const int iend_pack         = std::min( iend - istart,
                                                istart_pack + packsize );
        // const int current_pack_size = iend_pack - istart_pack;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                         \
                                      position_x /* [istart_pack:current_pack_size] */, \
                                      position_y /* [istart_pack:current_pack_size] */, \
                                      position_z /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  deltaold [0:3 * nparts], \
                                  Sx0 [0:kTmpArraySize],   \
                                  Sy0 [0:kTmpArraySize],   \
                                  Sz0 [0:kTmpArraySize],   \
                                  Sx1 [0:kTmpArraySize],   \
                                  Sy1 [0:kTmpArraySize],   \
                                  Sz1 [0:kTmpArraySize],   \
                                  DSx [0:kTmpArraySize],   \
                                  DSy [0:kTmpArraySize],   \
                                  DSz [0:kTmpArraySize] )  \
        deviceptr( position_x,                             \
                   position_y,                             \
                   position_z )

    // #pragma acc parallel present( iold [0:3 * nparts],      \
    //                               deltaold [0:3 * nparts] ) \
    //     deviceptr( position_x,                              \
    //                position_y,                              \
    //                position_z,                              \
    //                Sx0,                                     \
    //                Sy0,                                     \
    //                Sz0,                                     \
    //                DSx,                                     \
    //                DSy,                                     \
    //                DSz )

    #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            // Offset ipart to map [0, current_pack_size)
            const int ipart_pack = ipart - istart_pack;

            // -------------------------------------
            // Variable declaration & initialization
            // -------------------------------------

            // arrays used for the Esirkepov projection method
            for( unsigned int i=0; i<5; i++ ) {
                Sx1[ipart_pack+i*packsize] = 0.;
                Sy1[ipart_pack+i*packsize] = 0.;
                Sz1[ipart_pack+i*packsize] = 0.;
            }

            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------

            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            double delta = deltaold[0*packsize+ipart];
            double delta2 = delta*delta;
            Sx0[ipart_pack+0*packsize] = 0.;
            Sx0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sx0[ipart_pack+2*packsize] = 0.75-delta2;
            Sx0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sx0[ipart_pack+4*packsize] = 0.;

            delta = deltaold[1*packsize+ipart];
            delta2 = delta*delta;
            Sy0[ipart_pack+0*packsize] = 0.;
            Sy0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sy0[ipart_pack+2*packsize] = 0.75-delta2;
            Sy0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sy0[ipart_pack+4*packsize] = 0.;

            delta = deltaold[2*packsize+ipart];
            delta2 = delta*delta;
            Sz0[ipart_pack+0*packsize] = 0.;
            Sz0[ipart_pack+1*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sz0[ipart_pack+2*packsize] = 0.75-delta2;
            Sz0[ipart_pack+3*packsize] = 0.5 * ( delta2+delta+0.25 );
            Sz0[ipart_pack+4*packsize] = 0.;

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            const double xpn = position_x[ ipart ] * dx_inv;
            const int ip = std::round( xpn );
            const int ipo = iold[0*packsize+ipart];
            const int ip_m_ipo = ip-ipo-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[ipart_pack+(ip_m_ipo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sx1[ipart_pack+(ip_m_ipo+2)*packsize] = 0.75-delta2;
            Sx1[ipart_pack+(ip_m_ipo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            const double ypn = position_y[ ipart ] * dy_inv;
            const int jp = std::round( ypn );
            const int jpo = iold[1*packsize+ipart];
            const int jp_m_jpo = jp-jpo-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[ipart_pack+(jp_m_jpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sy1[ipart_pack+(jp_m_jpo+2)*packsize] = 0.75-delta2;
            Sy1[ipart_pack+(jp_m_jpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            const double zpn = position_z[ ipart ] * dz_inv;
            const int kp = std::round( zpn );
            const int kpo = iold[2*packsize+ipart];
            const int kp_m_kpo = kp-kpo-k_domain_begin;
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            Sz1[ipart_pack+(kp_m_kpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sz1[ipart_pack+(kp_m_kpo+2)*packsize] = 0.75-delta2;
            Sz1[ipart_pack+(kp_m_kpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            // computes Esirkepov coefficients
            for( int i=0; i < 5; i++ ) {
                DSx[ipart_pack+i*packsize] = Sx1[ipart_pack+i*packsize] - Sx0[ipart_pack+i*packsize];
                DSy[ipart_pack+i*packsize] = Sy1[ipart_pack+i*packsize] - Sy0[ipart_pack+i*packsize];
                DSz[ipart_pack+i*packsize] = Sz1[ipart_pack+i*packsize] - Sz0[ipart_pack+i*packsize];
            }

            // ---------------------------
            // Calculate the total current
            // ---------------------------

            iold[ipart+0*packsize] -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
            // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
            iold[ipart+1*packsize] -= 2;
            iold[ipart+2*packsize] -= 2;
        }

        // Jx^(d,p,p)
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( DSx [0:kTmpArraySize], sumX [0:kTmpArraySize] )

    // #pragma acc parallel deviceptr( DSx, sumX )

    #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            sumX[ipart_pack+0*packsize] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSx[ ipart_pack+(k-1)*packsize ];
            }
        }

        const int    z_size0                  = nprimz;
        const int    yz_size0                 = nprimz * nprimy;
        const double dx_ov_dt_inv_cell_volume = dx_ov_dt * inv_cell_volume;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                     \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  Jx [0:Jx_size],         \
                                  Sy0 [0:kTmpArraySize],   \
                                  Sz0 [0:kTmpArraySize],   \
                                  DSy [0:kTmpArraySize],   \
                                  DSz [0:kTmpArraySize],   \
                                  sumX [0:kTmpArraySize] ) \
        deviceptr( charge, weight ) vector_length( 8 )

    // #pragma acc parallel present( iold [0:3 * nparts], \
    //                               Jx [0:Jx_size] )    \
    //     deviceptr( charge, weight, Sy0,                \
    //                Sz0, DSy, DSz, sumX ) vector_length( 8 )

    #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            const double crx_p = dx_ov_dt_inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            const int linindex0 = iold[ipart+0*packsize]*yz_size0+iold[ipart+1*packsize]*z_size0+iold[ipart+2*packsize];
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( int k=0 ; k<5 ; k++ ) {
                for( int j=0 ; j<5 ; j++ ) {
                    const double tmp = crx_p * ( Sy0[ipart_pack+j*packsize]*Sz0[ipart_pack+k*packsize] +
                                                 0.5*DSy[ipart_pack+j*packsize]*Sz0[ipart_pack+k*packsize] +
                                                 0.5*DSz[ipart_pack+k*packsize]*Sy0[ipart_pack+j*packsize] +
                                                 one_third*DSy[ipart_pack+j*packsize]*DSz[ipart_pack+k*packsize] );
                    const int idx = linindex0 + j*z_size0 + k;
                    for( int i=1 ; i<5 ; i++ ) {
                        const double val = sumX[ipart_pack+(i)*packsize] * tmp;
                        const int    jdx = idx + i * yz_size0;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jx [ jdx ] += val;
                    }
                }
            }//i
        }

        // Jy^(p,d,p)
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( DSy [0:kTmpArraySize], \
                                  sumX [0:kTmpArraySize] )

    // #pragma acc parallel deviceptr( DSy, sumX )

    #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            sumX[ipart_pack+0*packsize] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSy[ ipart_pack+(k-1)*packsize ];
            }
        }

        const int    z_size1                  = nprimz;
        const int    yz_size1                 = nprimz * ( nprimy + 1 );
        const double dy_ov_dt_inv_cell_volume = dy_ov_dt * inv_cell_volume;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                     \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  Jy [0:Jy_size],         \
                                  Sx0 [0:kTmpArraySize],   \
                                  Sz0 [0:kTmpArraySize],   \
                                  DSx [0:kTmpArraySize],   \
                                  DSz [0:kTmpArraySize],   \
                                  sumX [0:kTmpArraySize] ) \
        deviceptr( charge, weight ) vector_length( 8 )

    // #pragma acc parallel present( iold [0:3 * nparts], \
    //                               Jy [0:Jy_size] )    \
    //     deviceptr( charge, weight, Sx0,                \
    //                Sz0, DSx, DSz, sumX ) vector_length( 8 )

    #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            const double cry_p = dy_ov_dt_inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            const int linindex1 = iold[ipart+0*packsize]*yz_size1+iold[ipart+1*packsize]*z_size1+iold[ipart+2*packsize];
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( int k=0 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    const double tmp = cry_p * ( Sz0[ipart_pack+k*packsize]*Sx0[ipart_pack+i*packsize] +
                                                 0.5*DSz[ipart_pack+k*packsize]*Sx0[ipart_pack+i*packsize] +
                                                 0.5*DSx[ipart_pack+i*packsize]*Sz0[ipart_pack+k*packsize] +
                                                 one_third*DSz[ipart_pack+k*packsize]*DSx[ipart_pack+i*packsize] );
                    const int idx = linindex1 + i*yz_size1 + k;
                    for( int j=1 ; j<5 ; j++ ) {
                        const double val = sumX[ipart_pack+(j)*packsize] * tmp;
                        const int    jdx = idx + j * z_size1;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jy [ jdx ] += val;
                    }
                }
            }//i
        }

        // Jz^(p,p,d)
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( DSz [0:kTmpArraySize], \
                                  sumX [0:kTmpArraySize] )

    // #pragma acc parallel deviceptr( DSz, sumX )

    #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            sumX[ipart_pack+0*packsize] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[ipart_pack+k*packsize] = sumX[ipart_pack+(k-1)*packsize]-DSz[ ipart_pack+(k-1)*packsize ];
            }
        }

        const int    z_size2                  = nprimz + 1;
        const int    yz_size2                 = ( nprimz + 1 ) * nprimy;
        const double dz_ov_dt_inv_cell_volume = dz_ov_dt * inv_cell_volume;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                     \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  Jz [0:Jz_size],         \
                                  Sx0 [0:kTmpArraySize],   \
                                  Sy0 [0:kTmpArraySize],   \
                                  DSx [0:kTmpArraySize],   \
                                  DSy [0:kTmpArraySize],   \
                                  sumX [0:kTmpArraySize] ) \
        deviceptr( charge, weight )

    // #pragma acc parallel present( iold [0:3 * nparts], \
    //                               Jz [0:Jz_size] )    \
    //     deviceptr( charge, weight, Sx0,                \
    //                Sy0, DSx, DSy, sumX ) vector_length( 8 )

    #pragma acc loop gang worker
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            const int ipart_pack = ipart - istart_pack;

            const double crz_p = dz_ov_dt_inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            const int linindex2 = iold[ipart+0*packsize]*yz_size2+iold[ipart+1*packsize]*z_size2+iold[ipart+2*packsize];
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( int k=1 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    for( int j=0 ; j<5 ; j++ ) {
                        const double tmp = crz_p * ( Sx0[ipart_pack+i*packsize]*Sy0[ipart_pack+j*packsize] +
                                                     0.5*DSx[ipart_pack+i*packsize]*Sy0[ipart_pack+j*packsize] +
                                                     0.5*DSy[ipart_pack+j*packsize]*Sx0[ipart_pack+i*packsize] +
                                                     one_third*DSx[ipart_pack+i*packsize]*DSy[ipart_pack+j*packsize] );
                        const int idx = linindex2 + j*z_size2 + i*yz_size2;
                        const double val = sumX[ipart_pack+(k)*packsize] * tmp;
                        const int    jdx = idx + k;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jz[ jdx ] += val;
                    }
                }
            }//i

        } // End for ipart

    } // End for ipack

    } // end currentDepositionKernel

    static inline void
    densityDepositionKernel3D( 
                                double *__restrict__ rho,
                                int rho_size,
                                const double *__restrict__ device_particle_position_x,
                                const double *__restrict__ device_particle_position_y,
                                const double *__restrict__ device_particle_position_z,
                                const short *__restrict__ device_particle_charge,
                                const double *__restrict__ device_particle_weight,
                                const int *__restrict__ host_bin_index,
                                unsigned int,
                                unsigned int,
                                unsigned int,
                                const double *__restrict__ invgf,
                                int *__restrict__ iold,
                                const double *__restrict__ deltaold,
                                const unsigned int    number_of_particles,
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

    const unsigned int bin_count      = 1;
    const int          nparts = host_bin_index[bin_count - 1];

    // TODO(Etienne M): Implement a cuda/hip kernel and enable particle 3D sorting/binning

    const double *const __restrict__ position_x = device_particle_position_x;
    const double *const __restrict__ position_y = device_particle_position_y;
    const double *const __restrict__ position_z = device_particle_position_z;
    const short  *const __restrict__ charge     = device_particle_charge;
    const double *const __restrict__ weight     = device_particle_weight;

    const double one_third = 1./3.;

    const istart = 0; 
    const iend   = nparts; 

    const int packsize = nparts;
    const int npack    = ( ( iend - istart ) + ( packsize - 1 ) ) / packsize; // divide + ceil npack.

    static constexpr bool kAutoDeviceFree = true;
    const std::size_t     kTmpArraySize   = 5 * packsize;

    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sx1{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sy1{ kTmpArraySize };
    smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> host_device_Sz1{ kTmpArraySize };

    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sx1 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sy1 );
    smilei::tools::gpu::HostDeviceMemoryManagement::DeviceAllocate( host_device_Sz1 );

    double *const __restrict__ Sx1  = host_device_Sx1.data();
    double *const __restrict__ Sy1  = host_device_Sy1.data();
    double *const __restrict__ Sz1  = host_device_Sz1.data();

    for (int ipack=0 ; ipack<npack ; ipack++) {
        const int istart_pack       = istart + ipack * packsize;
        const int iend_pack         = std::min( iend - istart,
                                                istart_pack + packsize );
        // const int current_pack_size = iend_pack - istart_pack;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                         \
                                      position_x /* [istart_pack:current_pack_size] */, \
                                      position_y /* [istart_pack:current_pack_size] */, \
                                      position_z /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  deltaold [0:3 * nparts], \
                                  Sx1 [0:kTmpArraySize],   \
                                  Sy1 [0:kTmpArraySize],   \
                                  Sz1 [0:kTmpArraySize] )  \
        deviceptr( position_x,                             \
                   position_y,                             \
                   position_z )

    // #pragma acc parallel present( iold [0:3 * nparts],      \
    //                               deltaold [0:3 * nparts] ) \
    //     deviceptr( position_x,                              \
    //                position_y,                              \
    //                position_z,                              \
    //                Sx0,                                     \
    //                Sy0,                                     \
    //                Sz0,                                     \
    //                DSx,                                     \
    //                DSy,                                     \
    //                DSz )

    #pragma acc loop gang worker vector
#endif
        for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
            // Offset ipart to map [0, current_pack_size)
            const int ipart_pack = ipart - istart_pack;

            // -------------------------------------
            // Variable declaration & initialization
            // -------------------------------------

            // arrays used for the Esirkepov projection method
            for( unsigned int i=0; i<5; i++ ) {
                Sx1[ipart_pack+i*packsize] = 0.;
                Sy1[ipart_pack+i*packsize] = 0.;
                Sz1[ipart_pack+i*packsize] = 0.;
            }

            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------

            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            double delta  ;
            double delta2 ;

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            const double xpn = position_x[ ipart ] * dx_inv;
            const int ip = std::round( xpn );
            const int ipo = iold[0*packsize+ipart];
            const int ip_m_ipo = ip-ipo-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[ipart_pack+(ip_m_ipo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sx1[ipart_pack+(ip_m_ipo+2)*packsize] = 0.75-delta2;
            Sx1[ipart_pack+(ip_m_ipo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            const double ypn = position_y[ ipart ] * dy_inv;
            const int jp = std::round( ypn );
            const int jpo = iold[1*packsize+ipart];
            const int jp_m_jpo = jp-jpo-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[ipart_pack+(jp_m_jpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sy1[ipart_pack+(jp_m_jpo+2)*packsize] = 0.75-delta2;
            Sy1[ipart_pack+(jp_m_jpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            const double zpn = position_z[ ipart ] * dz_inv;
            const int kp = std::round( zpn );
            const int kpo = iold[2*packsize+ipart];
            const int kp_m_kpo = kp-kpo-k_domain_begin;
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            Sz1[ipart_pack+(kp_m_kpo+1)*packsize] = 0.5 * ( delta2-delta+0.25 );
            Sz1[ipart_pack+(kp_m_kpo+2)*packsize] = 0.75-delta2;
            Sz1[ipart_pack+(kp_m_kpo+3)*packsize] = 0.5 * ( delta2+delta+0.25 );

            // ---------------------------
            // Calculate the total current
            // ---------------------------

            // iold[ipart+0*packsize] -= 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
            // // i/j/kpo stored with - i/j/k_domain_begin in Interpolator
            // iold[ipart+1*packsize] -= 2;
            // iold[ipart+2*packsize] -= 2;
        }

        //if (diag_flag) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                     \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
          #pragma acc parallel present( iold [0:3 * nparts], \
                                  rho [0:rho_size],          \
                                  Sx1 [0:kTmpArraySize],     \
                                  Sy1 [0:kTmpArraySize],     \
                                  Sz1 [0:kTmpArraySize])     \
                                  deviceptr( charge, weight )
          #pragma acc loop gang worker
#endif
          for( int ipart=istart_pack ; ipart<iend_pack; ipart++ ) {
              int ipart_pack = ipart - ipack*packsize;

              double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
              int z_size2 =  nprimz;
              int yz_size2 =  nprimz*nprimy;
              int linindex2 = iold[ipart+0*nparts]*yz_size2+iold[ipart+1*nparts]*z_size2+iold[ipart+2*nparts];

              #pragma acc loop vector
              for( int k=1 ; k<5 ; k++ ) {
                   for( int i=0 ; i<5 ; i++ ) {
                       for( int j=0 ; j<5 ; j++ ) {
                           int idx = linindex2 + j*z_size2 + i*yz_size2;
                           int jdx = idx + k;

                           SMILEI_ACCELERATOR_ATOMIC
                           rho[ jdx ] += charge_weight * Sx1[ipart_pack+i*packsize]*Sy1[ipart_pack+j*packsize]*Sz1[ipart_pack+k*packsize];
                       }//j
                   }//i
               }//k
          //} // End for ipart
        } // if diag_flag

    } // End for ipack
    } // end currentDepositionKernel

} // namespace acc

#endif
