
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

    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )
        #include <cmath>

        #include "Tools.h"
        #include "gpu.h"
    // #elif defined( __CUDACC__ ) 
    //     #include "Params.h"
    //     #include "gpu.h"
    #elif defined( __HIP__ )
        #include <hip/hip_runtime.h>

        #include "Params.h"
        #include "gpu.h"
    #endif


    #if defined( PRIVATE_SMILEI_USE_OPENMP_PROJECTION_IMPLEMENTATION )

    #include "Projector3D2OrderGPUKernelAcc.h"

namespace naive {

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

    const int packsize = nparts;

        // const int current_pack_size = iend_pack - istart_pack;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                         \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */, \
                                      position_x /* [istart_pack:current_pack_size] */, \
                                      position_y /* [istart_pack:current_pack_size] */, \
                                      position_z /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_OPENACC_MODE )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  deltaold [0:3 * nparts], \
                                  Jx[0:Jx_size], \
                                  Jy[0:Jy_size], \
                                  Jz[0:Jz_size] \
                                  ) \
        deviceptr( position_x,                             \
                   position_y,                             \
                   position_z, charge, weight )

    // #pragma acc parallel present( iold [0:3 * nparts],      \
    //                               deltaold_ [0:3 * nparts] ) \
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
        for( int ipart=0 ; ipart<nparts; ipart++ ) {

            // -------------------------------------
            // Variable declaration & initialization
            // -------------------------------------

            double Sx0[5];
            double Sx1[5];
            double Sy0[5];
            double Sy1[5];
            double Sz0[5];
            double Sz1[5];

            double DSx[5];
            double DSy[5];
            double DSz[5];

            double sumX[5];

            // arrays used for the Esirkepov projection method
            for( unsigned int i=0; i<5; i++ ) {
                Sx1[i] = 0.;
                Sy1[i] = 0.;
                Sz1[i] = 0.;
            }

            // --------------------------------------------------------
            // Locate particles & Calculate Esirkepov coef. S, DS and W
            // --------------------------------------------------------

            // locate the particle on the primal grid at former time-step & calculate coeff. S0
            double delta = deltaold[0*packsize+ipart];
            double delta2 = delta*delta;
            Sx0[0] = 0.;
            Sx0[1] = 0.5 * ( delta2-delta+0.25 );
            Sx0[2] = 0.75-delta2;
            Sx0[3] = 0.5 * ( delta2+delta+0.25 );
            Sx0[4] = 0.;

            delta = deltaold[1*packsize+ipart];
            delta2 = delta*delta;
            Sy0[0] = 0.;
            Sy0[1] = 0.5 * ( delta2-delta+0.25 );
            Sy0[2] = 0.75-delta2;
            Sy0[3] = 0.5 * ( delta2+delta+0.25 );
            Sy0[4] = 0.;

            delta = deltaold[2*packsize+ipart];
            delta2 = delta*delta;
            Sz0[0] = 0.;
            Sz0[1] = 0.5 * ( delta2-delta+0.25 );
            Sz0[2] = 0.75-delta2;
            Sz0[3] = 0.5 * ( delta2+delta+0.25 );
            Sz0[4] = 0.;

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            const double xpn = position_x[ ipart ] * dx_inv;
            const int ip = std::round( xpn );
            const int ip_m_ipo = ip-iold[0*packsize+ipart]-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[(ip_m_ipo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sx1[(ip_m_ipo+2)] = 0.75-delta2;
            Sx1[(ip_m_ipo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double ypn = position_y[ ipart ] * dy_inv;
            const int jp = std::round( ypn );
            const int jp_m_jpo = jp-iold[1*packsize+ipart]-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[(jp_m_jpo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sy1[(jp_m_jpo+2)] = 0.75-delta2;
            Sy1[(jp_m_jpo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double zpn = position_z[ ipart ] * dz_inv;
            const int kp = std::round( zpn );
            const int kp_m_kpo = kp-iold[2*packsize+ipart]-k_domain_begin;
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            Sz1[(kp_m_kpo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sz1[(kp_m_kpo+2)] = 0.75-delta2;
            Sz1[(kp_m_kpo+3)] = 0.5 * ( delta2+delta+0.25 );

            // computes Esirkepov coefficients
            for( int i=0; i < 5; i++ ) {
                DSx[i] = Sx1[i] - Sx0[i];
                DSy[i] = Sy1[i] - Sy0[i];
                DSz[i] = Sz1[i] - Sz0[i];
            }

            // ---------------------------
            // Calculate the total current
            // ---------------------------

            int ipo = iold[ipart+0*nparts] - 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
            // i/j/kpo stored with - i/j/k_domain_begin_ in Interpolator
            int jpo = iold[ipart+1*nparts] - 2;
            int kpo = iold[ipart+2*nparts] - 2;

            const int    z_size0                  = nprimz;
            const int    yz_size0                 = nprimz * nprimy;

            const int    z_size1                  = nprimz;
            const int    yz_size1                 = nprimz * ( nprimy + 1 );

            const int    z_size2                  = nprimz + 1;
            const int    yz_size2                 = ( nprimz + 1 ) * nprimy;

            // Jx^(d,p,p)

            sumX[0] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[k] = sumX[(k-1)]-DSx[ (k-1) ];
            }
                    
            const double crx_p = dx_ov_dt * inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            const int linindex0 = ipo*yz_size0 + jpo*nprimz + kpo;

            for( int k=0 ; k<5 ; k++ ) {
                for( int j=0 ; j<5 ; j++ ) {
                    const double tmp = crx_p * ( Sy0[j]*Sz0[k] +
                                                 0.5*DSy[j]*Sz0[k] +
                                                 0.5*DSz[k]*Sy0[j] +
                                                 one_third*DSy[j]*DSz[k] );
                    const int idx = linindex0 + j*nprimz + k;
                    for( int i=1 ; i<5 ; i++ ) {
                        const double val = sumX[i] * tmp;
                        const int    jdx = idx + i * yz_size0;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jx [ jdx ] += val;
                    }
                }
            }//i

            // Jy^(p,d,p)
            sumX[0] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[k] = sumX[k-1]-DSy[ k-1 ];
            }

            const double cry_p = dy_ov_dt * inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            //const int linindex1 = iold[0]* (nprimz * ( nprimy + 1 )) +iold[1]*nprimz+iold[2];
            const int linindex1 = ipo * yz_size1+ jpo * z_size1 + kpo;

            for( int k=0 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    const double tmp = cry_p * ( Sz0[k]*Sx0[i] +
                                                 0.5*DSz[k]*Sx0[i] +
                                                 0.5*DSx[i]*Sz0[k] +
                                                 one_third*DSz[k]*DSx[i] );
                    const int idx = linindex1 + i* (nprimz * ( nprimy + 1 )) + k;
                    for( int j=1 ; j<5 ; j++ ) {
                        const double val = sumX[j] * tmp;
                        const int    jdx = idx + j * nprimz;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jy [ jdx ] += val;
                    }
                }
            }


            // Jz^(p,p,d)

            sumX[0] = 0.;
            for( int k=1 ; k<5 ; k++ ) {
                sumX[k] = sumX[(k-1)]-DSz[ (k-1) ];
            }

            const double crz_p = dz_ov_dt * inv_cell_volume * static_cast<double>( charge[ipart] ) * weight[ipart];

            const int linindex2 = ipo * yz_size2 + jpo * z_size2 + kpo;
//            const int linindex2 = iold[0]*yz_size2+iold[1]*z_size2+iold[2];

            for( int k=1 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    for( int j=0 ; j<5 ; j++ ) {
                        const double tmp = crz_p * ( Sx0[i]*Sy0[j] +
                                                     0.5*DSx[i]*Sy0[j] +
                                                     0.5*DSy[j]*Sx0[i] +
                                                     one_third*DSx[i]*DSy[j] );
                        const int idx = linindex2 + j*z_size2 + i*yz_size2;
                        const double val = sumX[(k)] * tmp;
                        const int    jdx = idx + k;

                        SMILEI_ACCELERATOR_ATOMIC
                        Jz[ jdx ] += val;
                    }
                }
            }

        }

    } // end currentDepositionKernel

    static inline void
    currentAndDensityDepositionKernel3D( double *__restrict__ Jx,
                                         double *__restrict__ Jy,
                                         double *__restrict__ Jz,
                                         double *__restrict__ rho,
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
#elif defined( SMILEI_OPENACC_MODE )
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
#elif defined( SMILEI_OPENACC_MODE )
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
#elif defined( SMILEI_OPENACC_MODE )
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
#ifdef SMILEI_OPENACC_MODE
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
#elif defined( SMILEI_OPENACC_MODE )
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
#elif defined( SMILEI_OPENACC_MODE )
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
#ifdef SMILEI_OPENACC_MODE
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
#elif defined( SMILEI_OPENACC_MODE )
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
#elif defined( SMILEI_OPENACC_MODE )
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
#ifdef SMILEI_OPENACC_MODE
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

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp atomic update
#elif defined( SMILEI_OPENACC_MODE )
    #pragma acc atomic
#endif
                        Jz[ jdx ] += val;
                    }
                }
            }//i

        } // End for ipart

        //if (diag_flag) {
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                     \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_OPENACC_MODE )
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
#if defined( SMILEI_ACCELERATOR_GPU_OMP )
                           #pragma omp atomic update
#elif defined( SMILEI_OPENACC_MODE )
                           #pragma acc atomic
#endif
                           rho[ jdx ] += charge_weight * Sx1[ipart_pack+i*packsize]*Sy1[ipart_pack+j*packsize]*Sz1[ipart_pack+k*packsize];
                       }//j
                   }//i
               }//k
          //} // End for ipart
        } // if diag_flag

    } // End for ipack
    } // end currentDepositionKernel


} // namespace naive

    #else

namespace hip {
    namespace detail {
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

        template <typename ComputeFloat,
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize>
        __global__ void
        // __launch_bounds__(kWorkgroupSize, 1)
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
                                         ComputeFloat dx_inv,
                                         ComputeFloat dy_inv,
                                         ComputeFloat dz_inv,
                                         ComputeFloat dx_ov_dt,
                                         ComputeFloat dy_ov_dt,
                                         ComputeFloat dz_ov_dt,
                                         int          i_domain_begin,
                                         int          j_domain_begin,
                                         int          k_domain_begin,
                                         int          nprimy,
                                         int          nprimz,
                                         int          not_spectral )
        {
            // TODO(Etienne M): refactor this function. Break it into smaller
            // pieces (lds init/store, coeff computation, deposition etc..)
            // TODO(Etienne M): __ldg could be used to slightly improve GDS load
            // speed. This would only have an effect on Nvidia cards as this
            // operation is a no op on AMD.
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x * gridDim.y * gridDim.z;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int z_cluster_coordinate          = blockIdx.z;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y * gridDim.z + y_cluster_coordinate * gridDim.z + z_cluster_coordinate; // The indexing order is: x * ywidth * zwidth + y * zwidth + z
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );
            const unsigned int global_z_scratch_space_coordinate_offset = z_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );

            const int    GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 3 /* 3D */, 2 /* 2nd order interpolation */ );
            ComputeFloat one_third             = 1. / 3.;

            // NOTE: We gain from the particles not being sorted inside a
            // cluster because it reduces the bank conflicts one gets when
            // multiple threads access the same part of the shared memory. Such
            // "conflicted" accesses are serialized !
            // NOTE: We use a bit to much LDS. For Jx, the first row could be
            // discarded, for Jy we could remove the first column.

            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 3 /* 3D */, 2 /* 2nd order interpolation */ );

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
                const int *const __restrict__ iold        = &device_iold[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];
                ComputeFloat Sz0[5];
                ComputeFloat Sz1[5];
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
                {
                    const ComputeFloat delta  = deltaold[2 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sz0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sz0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sz0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz0[4] = static_cast<ComputeFloat>( 0.0 );
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
                {
                    // const int    kp             = static_cast<int>( zpn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat zpn      = static_cast<ComputeFloat>( device_particle_position_z[particle_index] ) * dz_inv;
                    const int          kp       = std::round( zpn );
                    const int          kpo      = iold[2 * particle_count];
                    const int          kp_m_kpo = kp - kpo - k_domain_begin;
                    const ComputeFloat delta    = zpn - static_cast<ComputeFloat>( kp );
                    const ComputeFloat delta2   = delta * delta;

                    Sz1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sz1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sz1[2] = 0.0; // Always set below
                    Sz1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sz1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sz1[kp_m_kpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz1[kp_m_kpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sz1[kp_m_kpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }

                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * dz_ov_dt;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int kpo = iold[2 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_z_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;

                // Jx

                for( unsigned int j = 0; j < 5; ++j ) {
                    for( unsigned int k = 0; k < 5; ++k ) {
                        ComputeFloat tmp = crx_p * (   Sy0[j]*Sz0[k] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sy1[j] - Sy0[j] )*Sz0[k] + ( Sz1[k] - Sz0[k] )*Sy0[j] ) 
                                                     +                           one_third  *   ( Sy1[j] - Sy0[j] )    *     ( Sz1[k] - Sz0[k] ) );
                        ComputeFloat tmp_reduction{};
                        const int jk_loc = ( ipo * GPUClusterWithGCWidth + jpo + j ) * GPUClusterWithGCWidth + kpo + k;
                        for( unsigned int i = 1; i < 5; ++i ) {
                            tmp_reduction -= ( Sx1[i-1] - Sx0[i-1] ) * tmp;
                            const int loc = i*GPUClusterWithGCWidth*GPUClusterWithGCWidth + jk_loc;
                            atomic::LDS::AddNoReturn( &Jx_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // Jy

                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int k = 0; k < 5; ++k ) {
                        ComputeFloat tmp = cry_p * (   Sx0[i]*Sz0[k] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sx1[i] - Sx0[i] )*Sz0[k] + ( Sz1[k] - Sz0[k] )*Sx0[i] ) 
                                                     +                           one_third  *   ( Sx1[i] - Sx0[i] )    *     ( Sz1[k] - Sz0[k] ) );
                        ComputeFloat tmp_reduction{};
                        const int ik_loc = (( i + ipo ) * GPUClusterWithGCWidth + jpo ) * GPUClusterWithGCWidth + kpo + k;
                        for( unsigned int j = 1; j < 5; ++j ) {
                            tmp_reduction -= ( Sy1[j-1] - Sy0[j-1] ) * tmp;
                            const int loc = j*GPUClusterWithGCWidth + ik_loc;
                            atomic::LDS::AddNoReturn( &Jy_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // Jz

                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int j = 0; j < 5; ++j ) {
                        ComputeFloat tmp = crz_p * (   Sx0[i]*Sy0[j] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sx1[i] - Sx0[i] )*Sy0[j] + ( Sy1[j] - Sy0[j] )*Sx0[i] ) 
                                                     +                           one_third  *   ( Sx1[i] - Sx0[i] )    *     ( Sy1[j] - Sy0[j] ) );
                        ComputeFloat tmp_reduction{};
                        const int ij_loc = (( i + ipo ) * GPUClusterWithGCWidth + (jpo + j)) * GPUClusterWithGCWidth + kpo;
                        for( unsigned int k = 1; k < 5; ++k ) {
                            tmp_reduction -= ( Sz1[k-1] - Sz0[k-1] ) * tmp;
                            const int loc =  k  + ij_loc;
                            atomic::LDS::AddNoReturn( &Jz_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }
            }

            __syncthreads();

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int local_x_scratch_space_coordinate = field_index / ( GPUClusterWithGCWidth * GPUClusterWithGCWidth );
                const unsigned int local_y_scratch_space_coordinate = ( field_index % ( GPUClusterWithGCWidth * GPUClusterWithGCWidth ) ) / GPUClusterWithGCWidth;
                const unsigned int local_z_scratch_space_coordinate = field_index % GPUClusterWithGCWidth;

                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate;
                const unsigned int global_y_scratch_space_coordinate = global_y_scratch_space_coordinate_offset + local_y_scratch_space_coordinate;
                const unsigned int global_z_scratch_space_coordinate = global_z_scratch_space_coordinate_offset + local_z_scratch_space_coordinate;

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int global_memory_index = ( global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate ) * nprimz + global_z_scratch_space_coordinate;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index],                                                                                             static_cast<double>( Jx_scratch_space[field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + /* We handle the FTDT/picsar */ not_spectral * global_x_scratch_space_coordinate * nprimz], static_cast<double>( Jy_scratch_space[field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index + /* We handle the FTDT/picsar */ not_spectral * (global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate)],                                                                                             static_cast<double>(  Jz_scratch_space[field_index] ) );
            }
        } // end DepositCurrent


        template <typename ComputeFloat,
                  typename ReductionFloat,
                  std::size_t kWorkgroupSize>
        __global__ void
        // __launch_bounds__(kWorkgroupSize, 1)
        DepositCurrentAndDensity_3D_Order2( double *__restrict__ device_Jx,
                                            double *__restrict__ device_Jy,
                                            double *__restrict__ device_Jz,
                                            double *__restrict__ device_rho,
                                            int Jx_size,
                                            int Jy_size,
                                            int Jz_size,
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
                                            ComputeFloat dx_inv,
                                            ComputeFloat dy_inv,
                                            ComputeFloat dz_inv,
                                            ComputeFloat dx_ov_dt,
                                            ComputeFloat dy_ov_dt,
                                            ComputeFloat dz_ov_dt,
                                            int          i_domain_begin,
                                            int          j_domain_begin,
                                            int          k_domain_begin,
                                            int          nprimy,
                                            int          nprimz,
                                            int          not_spectral )
        {
            // TODO(Etienne M): refactor this function. Break it into smaller
            // pieces (lds init/store, coeff computation, deposition etc..)
            // TODO(Etienne M): __ldg could be used to slightly improve GDS load
            // speed. This would only have an effect on Nvidia cards as this
            // operation is a no op on AMD.
            const unsigned int workgroup_size = kWorkgroupSize; // blockDim.x;
            const unsigned int bin_count      = gridDim.x * gridDim.y * gridDim.z;
            const unsigned int loop_stride    = workgroup_size; // This stride should enable better memory access coalescing

            const unsigned int x_cluster_coordinate          = blockIdx.x;
            const unsigned int y_cluster_coordinate          = blockIdx.y;
            const unsigned int z_cluster_coordinate          = blockIdx.z;
            const unsigned int workgroup_dedicated_bin_index = x_cluster_coordinate * gridDim.y * gridDim.z + y_cluster_coordinate * gridDim.z + z_cluster_coordinate; // The indexing order is: x * ywidth * zwidth + y * zwidth + z
            const unsigned int thread_index_offset           = threadIdx.x;

            // The unit is the cell
            const unsigned int global_x_scratch_space_coordinate_offset = x_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );
            const unsigned int global_y_scratch_space_coordinate_offset = y_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );
            const unsigned int global_z_scratch_space_coordinate_offset = z_cluster_coordinate * Params::getGPUClusterWidth( 3 /* 3D */ );

            const int    GPUClusterWithGCWidth = Params::getGPUClusterWithGhostCellWidth( 3 /* 3D */, 2 /* 2nd order interpolation */ );
            ComputeFloat one_third             = 1. / 3.;

            // NOTE: We gain from the particles not being sorted inside a
            // cluster because it reduces the bank conflicts one gets when
            // multiple threads access the same part of the shared memory. Such
            // "conflicted" accesses are serialized !
            // NOTE: We use a bit to much LDS. For Jx, the first row could be
            // discarded, for Jy we could remove the first column.

            static constexpr unsigned int kFieldScratchSpaceSize = Params::getGPUInterpolationClusterCellVolume( 3 /* 3D */, 2 /* 2nd order interpolation */ );

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
                const int *const __restrict__ iold        = &device_iold[particle_index];
                const double *const __restrict__ deltaold = &device_deltaold_[particle_index];

                ComputeFloat Sx0[5];
                ComputeFloat Sx1[5];
                ComputeFloat Sy0[5];
                ComputeFloat Sy1[5];
                ComputeFloat Sz0[5];
                ComputeFloat Sz1[5];
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
                {
                    const ComputeFloat delta  = deltaold[2 * particle_count];
                    const ComputeFloat delta2 = delta * delta;

                    Sz0[0] = static_cast<ComputeFloat>( 0.0 );
                    Sz0[1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz0[2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sz0[3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz0[4] = static_cast<ComputeFloat>( 0.0 );
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
                {
                    // const int    kp             = static_cast<int>( zpn + 0.5 ); // std::round | rounding approximation which is correct enough and faster in this case
                    const ComputeFloat zpn      = static_cast<ComputeFloat>( device_particle_position_z[particle_index] ) * dz_inv;
                    const int          kp       = std::round( zpn );
                    const int          kpo      = iold[2 * particle_count];
                    const int          kp_m_kpo = kp - kpo - k_domain_begin;
                    const ComputeFloat delta    = zpn - static_cast<ComputeFloat>( kp );
                    const ComputeFloat delta2   = delta * delta;

                    Sz1[0] = static_cast<ComputeFloat>( 0.0 );
                    Sz1[1] = static_cast<ComputeFloat>( 0.0 );
                    // Sz1[2] = 0.0; // Always set below
                    Sz1[3] = static_cast<ComputeFloat>( 0.0 );
                    Sz1[4] = static_cast<ComputeFloat>( 0.0 );

                    Sz1[kp_m_kpo + 1] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 - delta + static_cast<ComputeFloat>( 0.25 ) );
                    Sz1[kp_m_kpo + 2] = static_cast<ComputeFloat>( 0.75 ) - delta2;
                    Sz1[kp_m_kpo + 3] = static_cast<ComputeFloat>( 0.5 ) * ( delta2 + delta + static_cast<ComputeFloat>( 0.25 ) );
                }

                // (x,y,z) components of the current density for the macro-particle
                const ComputeFloat charge_weight = inv_cell_volume * static_cast<ComputeFloat>( device_particle_charge[particle_index] ) * static_cast<ComputeFloat>( device_particle_weight[particle_index] );
                const ComputeFloat crx_p         = charge_weight * dx_ov_dt;
                const ComputeFloat cry_p         = charge_weight * dy_ov_dt;
                const ComputeFloat crz_p         = charge_weight * dz_ov_dt;

                // This is the particle position as grid index
                // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
                const int ipo = iold[0 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_x_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int jpo = iold[1 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_y_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;
                const int kpo = iold[2 * particle_count] -
                                2 /* Offset so we dont uses negative numbers in the loop */ -
                                global_z_scratch_space_coordinate_offset /* Offset to get cluster relative coordinates */;

                // Jx

                for( unsigned int j = 0; j < 5; ++j ) {
                    for( unsigned int k = 0; k < 5; ++k ) {
                        ComputeFloat tmp = crx_p * (   Sy0[j]*Sz0[k] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sy1[j] - Sy0[j] )*Sz0[k] + ( Sz1[k] - Sz0[k] )*Sy0[j] ) 
                                                     +                           one_third  *   ( Sy1[j] - Sy0[j] )    *     ( Sz1[k] - Sz0[k] ) );
                        ComputeFloat tmp_reduction{};
                        const int jk_loc = ( ipo  * GPUClusterWithGCWidth + jpo + j ) * GPUClusterWithGCWidth + kpo + k;
                        for( unsigned int i = 1; i < 5; ++i ) {
                            tmp_reduction -= ( Sx1[i-1] - Sx0[i-1] ) * tmp;
                            const int loc = i*GPUClusterWithGCWidth*GPUClusterWithGCWidth + jk_loc;
                            atomic::LDS::AddNoReturn( &Jx_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // Jy

                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int k = 0; k < 5; ++k ) {
                        ComputeFloat tmp = cry_p * (   Sx0[i]*Sz0[k] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sx1[i] - Sx0[i] )*Sz0[k] + ( Sz1[k] - Sz0[k] )*Sx0[i] ) 
                                                     +                           one_third  *   ( Sx1[i] - Sx0[i] )    *     ( Sz1[k] - Sz0[k] ) );
                        ComputeFloat tmp_reduction{};
                        const int ik_loc = (( i + ipo ) * GPUClusterWithGCWidth + jpo ) * GPUClusterWithGCWidth + kpo + k;
                        for( unsigned int j = 1; j < 5; ++j ) {
                            tmp_reduction -= ( Sy1[j-1] - Sy0[j-1] ) * tmp;
                            const int loc = j*GPUClusterWithGCWidth + ik_loc;
                            atomic::LDS::AddNoReturn( &Jy_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }

                // Jz

                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int j = 0; j < 5; ++j ) {
                        ComputeFloat tmp = crz_p * (   Sx0[i]*Sy0[j] 
                                                     + static_cast<ComputeFloat>( 0.5 )     * ( ( Sx1[i] - Sx0[i] )*Sy0[j] + ( Sy1[j] - Sy0[j] )*Sx0[i] ) 
                                                     +                           one_third  *   ( Sx1[i] - Sx0[i] )    *     ( Sy1[j] - Sy0[j] ) );
                        ComputeFloat tmp_reduction{};
                        const int ij_loc = (( i + ipo ) * GPUClusterWithGCWidth + jpo + j) * GPUClusterWithGCWidth + kpo;
                        for( unsigned int k = 1; k < 5; ++k ) {
                            tmp_reduction -= ( Sz1[k-1] - Sz0[k-1] ) * tmp;
                            const int loc =  ij_loc + k;
                            atomic::LDS::AddNoReturn( &Jz_scratch_space[loc], static_cast<ReductionFloat>( tmp_reduction ) );
                        }
                    }
                }


                // Rho

                for( unsigned int i = 0; i < 5; ++i ) {
                    for( unsigned int j = 0; j < 5; ++j ) {
                        ComputeFloat tmp = charge_weight * Sx1[i]*Sy1[j]; 
                        const int ij_loc = (( i + ipo ) * GPUClusterWithGCWidth + (jpo + j)) * GPUClusterWithGCWidth + kpo;
                        for( unsigned int k = 0; k < 5; ++k ) {
                            const int loc =  ij_loc + k;
                            atomic::LDS::AddNoReturn( &rho_scratch_space[loc], static_cast<ReductionFloat>( tmp * Sz1[k] ) );
                        }
                    }
                }
            }

            __syncthreads();

            for( unsigned int field_index = thread_index_offset;
                 field_index < kFieldScratchSpaceSize;
                 field_index += workgroup_size ) {

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int local_x_scratch_space_coordinate = field_index / ( GPUClusterWithGCWidth * GPUClusterWithGCWidth );
                const unsigned int local_y_scratch_space_coordinate = ( field_index % ( GPUClusterWithGCWidth * GPUClusterWithGCWidth ) ) / GPUClusterWithGCWidth;
                const unsigned int local_z_scratch_space_coordinate = field_index % GPUClusterWithGCWidth;

                const unsigned int global_x_scratch_space_coordinate = global_x_scratch_space_coordinate_offset + local_x_scratch_space_coordinate;
                const unsigned int global_y_scratch_space_coordinate = global_y_scratch_space_coordinate_offset + local_y_scratch_space_coordinate;
                const unsigned int global_z_scratch_space_coordinate = global_z_scratch_space_coordinate_offset + local_z_scratch_space_coordinate;

                // The indexing order is: x * ywidth * zwidth + y * zwidth + z
                const unsigned int global_memory_index = ( global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate ) * nprimz + global_z_scratch_space_coordinate;

                // These atomics are basically free (very few of them).
                atomic::GDS::AddNoReturn( &device_Jx[global_memory_index], static_cast<double>( Jx_scratch_space[field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jy[global_memory_index + /* We handle the FTDT/picsar */ not_spectral * global_x_scratch_space_coordinate * nprimz], static_cast<double>( Jy_scratch_space[field_index] ) );
                atomic::GDS::AddNoReturn( &device_Jz[global_memory_index + /* We handle the FTDT/picsar */ not_spectral * ( global_x_scratch_space_coordinate * nprimy + global_y_scratch_space_coordinate )], static_cast<double>( Jz_scratch_space[field_index] ) );
                atomic::GDS::AddNoReturn( &device_rho[global_memory_index], static_cast<double>( rho_scratch_space[field_index] ) );
            }
        }
    } // namespace kernel


    static inline void
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
                               int    not_spectral )
    {
        SMILEI_ASSERT( Params::getGPUClusterWidth( 3 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        // NOTE:
        // This cluster is very strongly bound by atomic operations in LDS (shared memory)
        // TODO(Etienne M): Find a way to lessen the atomic usage

        const ::dim3 kGridDimension /* In blocks */ { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), static_cast<uint32_t>( z_dimension_bin_count ) };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // NOTE: On cards lacking hardware backed Binary64 atomic operations,
        // falling back to Binary32 (supposing hardware support for atomic
        // operations) can lead to drastic performance improvement.
        // One just need to assign 'float' to ReductionFloat.
        //
        using ComputeFloat   = double;
        using ReductionFloat = double;

        auto KernelFunction = kernel::DepositCurrentDensity_3D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;

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
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv, dy_inv, dz_inv,
                            dx_ov_dt, dy_ov_dt, dz_ov_dt,
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral );

        checkHIPErrors( ::hipDeviceSynchronize() );
    }

    static inline void
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
        SMILEI_ASSERT( Params::getGPUClusterWidth( 3 /* 2D */ ) != -1 &&
                       Params::getGPUClusterGhostCellBorderWidth( 2 /* 2nd order interpolation */ ) != -1 );

        // NOTE:
        // This cluster is very strongly bound by atomic operations in LDS (shared memory)
        // TODO(Etienne M): Find a way to lessen the atomic usage

        const ::dim3 kGridDimension /* In blocks */ { static_cast<uint32_t>( x_dimension_bin_count ), static_cast<uint32_t>( y_dimension_bin_count ), static_cast<uint32_t>( z_dimension_bin_count ) };

        static constexpr std::size_t kWorkgroupSize = 128;
        const ::dim3                 kBlockDimension{ static_cast<uint32_t>( kWorkgroupSize ), 1, 1 };

        // NOTE: On cards lacking hardware backed Binary64 atomic operations,
        // falling back to Binary32 (supposing hardware support for atomic
        // operations) can lead to drastic performance improvement.
        // One just need to assign 'float' to ReductionFloat.
        //
        using ComputeFloat   = double;
        using ReductionFloat = double;

        auto KernelFunction = kernel::DepositCurrentAndDensity_3D_Order2<ComputeFloat, ReductionFloat, kWorkgroupSize>;

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
                            device_particle_position_z,
                            device_particle_charge,
                            device_particle_weight,
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_bin_index ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_invgf_ ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_iold ),
                            smilei::tools::gpu::HostDeviceMemoryManagement::GetDevicePointer( host_deltaold_ ),
                            inv_cell_volume,
                            dx_inv, dy_inv, dz_inv,
                            dx_ov_dt, dy_ov_dt, dz_ov_dt,
                            i_domain_begin, j_domain_begin, k_domain_begin,
                            nprimy, nprimz,
                            not_spectral );

        checkHIPErrors( ::hipDeviceSynchronize() );
    }

} // namespace hip

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
