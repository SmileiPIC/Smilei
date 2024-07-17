//! Naive ACC/OMP implementation

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

        // const int current_pack_size = iend_pack - istart_pack;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                         \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */, \
                                      position_x /* [istart_pack:current_pack_size] */, \
                                      position_y /* [istart_pack:current_pack_size] */, \
                                      position_z /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
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
            double delta = deltaold[0*nparts+ipart];
            double delta2 = delta*delta;
            Sx0[0] = 0.;
            Sx0[1] = 0.5 * ( delta2-delta+0.25 );
            Sx0[2] = 0.75-delta2;
            Sx0[3] = 0.5 * ( delta2+delta+0.25 );
            Sx0[4] = 0.;

            delta = deltaold[1*nparts+ipart];
            delta2 = delta*delta;
            Sy0[0] = 0.;
            Sy0[1] = 0.5 * ( delta2-delta+0.25 );
            Sy0[2] = 0.75-delta2;
            Sy0[3] = 0.5 * ( delta2+delta+0.25 );
            Sy0[4] = 0.;

            delta = deltaold[2*nparts+ipart];
            delta2 = delta*delta;
            Sz0[0] = 0.;
            Sz0[1] = 0.5 * ( delta2-delta+0.25 );
            Sz0[2] = 0.75-delta2;
            Sz0[3] = 0.5 * ( delta2+delta+0.25 );
            Sz0[4] = 0.;

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            const double xpn = position_x[ ipart ] * dx_inv;
            const int ip = std::round( xpn );
            const int ip_m_ipo = ip-iold[0*nparts+ipart]-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[(ip_m_ipo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sx1[(ip_m_ipo+2)] = 0.75-delta2;
            Sx1[(ip_m_ipo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double ypn = position_y[ ipart ] * dy_inv;
            const int jp = std::round( ypn );
            const int jp_m_jpo = jp-iold[1*nparts+ipart]-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[(jp_m_jpo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sy1[(jp_m_jpo+2)] = 0.75-delta2;
            Sy1[(jp_m_jpo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double zpn = position_z[ ipart ] * dz_inv;
            const int kp = std::round( zpn );
            const int kp_m_kpo = kp-iold[2*nparts+ipart]-k_domain_begin;
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

        } // end particle loop

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

        // const int current_pack_size = iend_pack - istart_pack;

#if defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target is_device_ptr( /* to: */                                         \
                                      charge /* [istart_pack:current_pack_size] */, \
                                      weight /* [istart_pack:current_pack_size] */, \
                                      position_x /* [istart_pack:current_pack_size] */, \
                                      position_y /* [istart_pack:current_pack_size] */, \
                                      position_z /* [istart_pack:current_pack_size] */ )
    #pragma omp teams distribute parallel for
#elif defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( iold [0:3 * nparts],     \
                                  deltaold [0:3 * nparts], \
                                  rho[0:rho_size] \
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

            double Sx1[5];
            double Sy1[5];
            double Sz1[5];

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
            double delta;
            double delta2;

            // locate the particle on the primal grid at current time-step & calculate coeff. S1
            const double xpn = position_x[ ipart ] * dx_inv;
            const int ip = std::round( xpn );
            const int ip_m_ipo = ip-iold[0*nparts+ipart]-i_domain_begin;
            delta  = xpn - ( double )ip;
            delta2 = delta*delta;
            Sx1[(ip_m_ipo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sx1[(ip_m_ipo+2)] = 0.75-delta2;
            Sx1[(ip_m_ipo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double ypn = position_y[ ipart ] * dy_inv;
            const int jp = std::round( ypn );
            const int jp_m_jpo = jp-iold[1*nparts+ipart]-j_domain_begin;
            delta  = ypn - ( double )jp;
            delta2 = delta*delta;
            Sy1[(jp_m_jpo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sy1[(jp_m_jpo+2)] = 0.75-delta2;
            Sy1[(jp_m_jpo+3)] = 0.5 * ( delta2+delta+0.25 );

            const double zpn = position_z[ ipart ] * dz_inv;
            const int kp = std::round( zpn );
            const int kp_m_kpo = kp-iold[2*nparts+ipart]-k_domain_begin;
            delta  = zpn - ( double )kp;
            delta2 = delta*delta;
            Sz1[(kp_m_kpo+1)] = 0.5 * ( delta2-delta+0.25 );
            Sz1[(kp_m_kpo+2)] = 0.75-delta2;
            Sz1[(kp_m_kpo+3)] = 0.5 * ( delta2+delta+0.25 );

            // ---------------------------
            // Calculate the total current
            // ---------------------------

            int ipo = iold[ipart+0*nparts] - 2;   //This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
            // i/j/kpo stored with - i/j/k_domain_begin_ in Interpolator
            int jpo = iold[ipart+1*nparts] - 2;
            int kpo = iold[ipart+2*nparts] - 2;

            // Rho

            double charge_weight = inv_cell_volume * ( double )( charge[ ipart ] )*weight[ ipart ];
            int z_size2 =  nprimz;
            int yz_size2 =  nprimz*nprimy;
            int linindex2 = ipo*yz_size2+jpo*z_size2+kpo;

            for( int k=1 ; k<5 ; k++ ) {
                for( int i=0 ; i<5 ; i++ ) {
                    for( int j=0 ; j<5 ; j++ ) {
                        int idx = linindex2 + j*z_size2 + i*yz_size2;
                        int jdx = idx + k;
                        SMILEI_ACCELERATOR_ATOMIC
                        rho[ jdx ] += charge_weight * Sx1[i]*Sy1[j]*Sz1[k];
                    } //j
                } //i
            } //k

        } // End for ipart

    } // end densityDepositionKernel

} // namespace acc

#endif
