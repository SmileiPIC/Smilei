#if defined( SMILEI_ACCELERATOR_GPU_OMP )

    #include <cmath>

    #include "Tools.h"

static inline void
naiveCurrentDepositionKernel( double *__restrict__ Jx,
                              double *__restrict__ Jy,
                              double *__restrict__ Jz,
                              int Jx_size,
                              int Jy_size,
                              int Jz_size,
                              const double *__restrict__ particle_position_x,
                              const double *__restrict__ particle_position_y,
                              const double *__restrict__ particle_momentum_z,
                              const short *__restrict__ particle_charge,
                              const double *__restrict__ particle_weight,
                              int *__restrict__ bin_index,
                              int bin_count,
                              const double *__restrict__ invgf_,
                              const int *__restrict__ iold_,
                              const double *__restrict__ deltaold_,
                              double inv_cell_volume,
                              double dx_inv,
                              double dy_inv,
                              double dx_ov_dt,
                              double dy_ov_dt,
                              int    i_domain_begin,
                              int    j_domain_begin,
                              int    nprimy,
                              int    pxr )
{
    SMILEI_ASSERT( bin_count > 0 );

    const int particle_count = bin_index[bin_count - 1];

    // // Arrays used for the Esirkepov projection method
    // static constexpr bool kAutoDeviceFree = true;
    // const std::size_t     kTmpArraySize   = particle_count * 5;

    // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> Sx0_buffer{ kTmpArraySize };
    // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> Sx1_buffer{ kTmpArraySize };
    // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> Sy0_buffer{ kTmpArraySize };
    // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> Sy1_buffer{ kTmpArraySize };
    // // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> DSx_buffer{ kTmpArraySize };
    // // smilei::tools::gpu::NonInitializingVector<double, kAutoDeviceFree> DSy_buffer{ kTmpArraySize };

    // double *const __restrict__ Sx0_buffer_data = Sx0_buffer.data();
    // double *const __restrict__ Sx1_buffer_data = Sx1_buffer.data();
    // double *const __restrict__ Sy0_buffer_data = Sy0_buffer.data();
    // double *const __restrict__ Sy1_buffer_data = Sy1_buffer.data();
    // // double *const __restrict__ DSx_buffer_data = DSx_buffer.data();
    // // double *const __restrict__ DSy_buffer_data = DSy_buffer.data();

    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
    const int interpolation_range_2D_size = particle_count + 1 * particle_count;

        #pragma omp target is_device_ptr /* map */ ( /* to: */                                     \
                                                     particle_position_x /* [0:particle_count] */, \
                                                     particle_position_y /* [0:particle_count] */, \
                                                     particle_momentum_z /* [0:particle_count] */, \
                                                     particle_charge /* [0:particle_count] */,     \
                                                     particle_weight /* [0:particle_count] */ )
        //        map( from                                            \
                    //  : Sx0_buffer_data [0:kTmpArraySize],            \ 
                    //    Sx1_buffer_data [0:kTmpArraySize],            \
                    //    Sy0_buffer_data [0:kTmpArraySize],            \
                    //    Sy1_buffer_data [0:kTmpArraySize] )
        #pragma omp teams      thread_limit( 64 )
        #pragma omp distribute parallel for
    #endif
    for( int particle_index = 0; particle_index < particle_count; ++particle_index ) {
        const double invgf                        = invgf_[particle_index];
        const int *const __restrict__ iold        = &iold_[particle_index];
        const double *const __restrict__ deltaold = &deltaold_[particle_index];

        double Sx0[5];
        double Sx1[5];
        double Sy0[5];
        double Sy1[5];
        // double DSx[5];
        // double DSy[5];

        // double *const __restrict__ Sx0 = Sx0_buffer_data + 5 * ( particle_index - 0 );
        // double *const __restrict__ Sx1 = Sx1_buffer_data + 5 * ( particle_index - 0 );
        // double *const __restrict__ Sy0 = Sy0_buffer_data + 5 * ( particle_index - 0 );
        // double *const __restrict__ Sy1 = Sy1_buffer_data + 5 * ( particle_index - 0 );
        // // double *const __restrict__ DSx = DSx_buffer_data + 5 * ( particle_index - 0 );
        // // double *const __restrict__ DSy = DSy_buffer_data + 5 * ( particle_index - 0 );

        // Variable declaration & initialization
        // Esirkepov's paper: https://arxiv.org/pdf/physics/9901047.pdf

        // Locate the particle on the primal grid at former time-step & calculate coeff. S0
        {
            const double delta  = deltaold[0 * particle_count];
            const double delta2 = delta * delta;
            Sx0[0]              = 0.0;
            Sx0[1]              = 0.5 * ( delta2 - delta + 0.25 );
            Sx0[2]              = 0.75 - delta2;
            Sx0[3]              = 0.5 * ( delta2 + delta + 0.25 );
            Sx0[4]              = 0.0;
        }
        {
            const double delta  = deltaold[1 * particle_count];
            const double delta2 = delta * delta;
            Sy0[0]              = 0.0;
            Sy0[1]              = 0.5 * ( delta2 - delta + 0.25 );
            Sy0[2]              = 0.75 - delta2;
            Sy0[3]              = 0.5 * ( delta2 + delta + 0.25 );
            Sy0[4]              = 0.0;
        }

        // Locate the particle on the primal grid at current time-step & calculate coeff. S1
        {
            const double xpn      = particle_position_x[particle_index] * dx_inv;
            const int    ip       = std::round( xpn );
            const int    ipo      = iold[0 * particle_count];
            const int    ip_m_ipo = ip - ipo - i_domain_begin;
            const double delta    = xpn - static_cast<double>( ip );
            const double delta2   = delta * delta;

            Sx1[0] = 0.0;
            Sx1[1] = 0.0;
            // Sx1[2] = 0.0; // Always set below
            Sx1[3] = 0.0;
            Sx1[4] = 0.0;

            Sx1[ip_m_ipo + 1] = 0.5 * ( delta2 - delta + 0.25 );
            Sx1[ip_m_ipo + 2] = 0.75 - delta2;
            Sx1[ip_m_ipo + 3] = 0.5 * ( delta2 + delta + 0.25 );
        }
        {
            const double ypn      = particle_position_y[particle_index] * dy_inv;
            const int    jp       = std::round( ypn );
            const int    jpo      = iold[1 * particle_count];
            const int    jp_m_jpo = jp - jpo - j_domain_begin;
            const double delta    = ypn - static_cast<double>( jp );
            const double delta2   = delta * delta;

            Sy1[0] = 0.0;
            Sy1[1] = 0.0;
            // Sy1[2] = 0.0; // Always set below
            Sy1[3] = 0.0;
            Sy1[4] = 0.0;

            Sy1[jp_m_jpo + 1] = 0.5 * ( delta2 - delta + 0.25 );
            Sy1[jp_m_jpo + 2] = 0.75 - delta2;
            Sy1[jp_m_jpo + 3] = 0.5 * ( delta2 + delta + 0.25 );
        }

        // DSx[0] = Sx1[0] - Sx0[0];
        // DSx[1] = Sx1[1] - Sx0[1];
        // DSx[2] = Sx1[2] - Sx0[2];
        // DSx[3] = Sx1[3] - Sx0[3];
        // DSx[4] = Sx1[4] - Sx0[4];

        // DSy[0] = Sy1[0] - Sy0[0];
        // DSy[1] = Sy1[1] - Sy0[1];
        // DSy[2] = Sy1[2] - Sy0[2];
        // DSy[3] = Sy1[3] - Sy0[3];
        // DSy[4] = Sy1[4] - Sy0[4];
        // }

        // // Charge deposition on the grid

        // for( int particle_index = 0; particle_index < particle_count; ++particle_index ) {
        //     const double invgf                        = invgf_[particle_index];
        //     const int *const __restrict__ iold        = &iold_[particle_index];
        //     const double *const __restrict__ deltaold = &deltaold_[particle_index];

        //     double *const __restrict__ Sx0 = Sx0_buffer_data + 5 * ( particle_index - 0 );
        //     double *const __restrict__ Sx1 = Sx1_buffer_data + 5 * ( particle_index - 0 );
        //     double *const __restrict__ Sy0 = Sy0_buffer_data + 5 * ( particle_index - 0 );
        //     double *const __restrict__ Sy1 = Sy1_buffer_data + 5 * ( particle_index - 0 );
        //     // double *const __restrict__ DSx = DSx_buffer_data + 5 * ( particle_index - 0 );
        //     // double *const __restrict__ DSy = DSy_buffer_data + 5 * ( particle_index - 0 );

        // (x,y,z) components of the current density for the macro-particle
        const double charge_weight = inv_cell_volume * static_cast<double>( particle_charge[particle_index] ) * particle_weight[particle_index];
        const double crx_p         = charge_weight * dx_ov_dt;
        const double cry_p         = charge_weight * dy_ov_dt;
        const double crz_p         = charge_weight * ( 1.0 / 3.0 ) * particle_momentum_z[particle_index] * invgf;

        // This is the particle position as grid index
        // This minus 2 come from the order 2 scheme, based on a 5 points stencil from -2 to +2.
        const int ipo = iold[0 * particle_count] - 2;
        const int jpo = iold[1 * particle_count] - 2;

        for( unsigned int i = 0; i < 1; ++i ) {
            const int iloc = ( i + ipo ) * nprimy + jpo;
            /* Jx[iloc] += tmpJx[0]; */
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
            Jz[iloc] += crz_p * ( Sy1[0] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
            double tmp = 0.0;
            for( unsigned int j = 1; j < 5; j++ ) {
                tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
                Jy[iloc + j + pxr * ( /* i + */ ipo )] += tmp;
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
                Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] /* + Sx0[i] */ ) +
                                          Sy1[j] * ( /* 0.5 * Sx0[i] + */ Sx1[i] ) );
            }
        }

        double tmpJx[5]{};

        for( unsigned int i = 1; i < 5; ++i ) {
            const int iloc = ( i + ipo ) * nprimy + jpo;
            tmpJx[0] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( 0.5 * ( Sy1[0] - Sy0[0] ) );
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
            Jx[iloc] += tmpJx[0];
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
            Jz[iloc] += crz_p * ( Sy1[0] * ( 0.5 * Sx0[i] + Sx1[i] ) );
            double tmp = 0.0;
            for( unsigned int j = 1; j < 5; ++j ) {
                tmpJx[j] -= crx_p * ( Sx1[i - 1] - Sx0[i - 1] ) * ( Sy0[j] + 0.5 * ( Sy1[j] - Sy0[j] ) );
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
                Jx[iloc + j] += tmpJx[j];
                tmp -= cry_p * ( Sy1[j - 1] - Sy0[j - 1] ) * ( Sx0[i] + 0.5 * ( Sx1[i] - Sx0[i] ) );
    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
                Jy[iloc + j + pxr * ( i + ipo )] += tmp;

    #if defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp atomic update
    #endif
                Jz[iloc + j] += crz_p * ( Sy0[j] * ( 0.5 * Sx1[i] + Sx0[i] ) +
                                          Sy1[j] * ( 0.5 * Sx0[i] + Sx1[i] ) );
            }
        }
    }
}

void hipCurrentDepositionKernel( double *__restrict__ Jx,
                                 double *__restrict__ Jy,
                                 double *__restrict__ Jz,
                                 int Jx_size,
                                 int Jy_size,
                                 int Jz_size,
                                 const double *__restrict__ particle_position_x,
                                 const double *__restrict__ particle_position_y,
                                 const double *__restrict__ particle_momentum_z,
                                 const short *__restrict__ particle_charge,
                                 const double *__restrict__ particle_weight,
                                 int *__restrict__ bin_index,
                                 int bin_count,
                                 const double *__restrict__ invgf_,
                                 const int *__restrict__ iold_,
                                 const double *__restrict__ deltaold_,
                                 double inv_cell_volume,
                                 double dx_inv,
                                 double dy_inv,
                                 double dx_ov_dt,
                                 double dy_ov_dt,
                                 int    i_domain_begin,
                                 int    j_domain_begin,
                                 int    nprimy,
                                 int    pxr )
{
}

// 3 streams (Jx Jy Jz)
//

// TODO(Etienne M): Change .cpp to .cu

//! Project global current densities (EMfields->Jx_/Jy_/Jz_)
//!
extern "C" void
currentDepositionKernel( double *__restrict__ Jx,
                         double *__restrict__ Jy,
                         double *__restrict__ Jz,
                         int Jx_size,
                         int Jy_size,
                         int Jz_size,
                         const double *__restrict__ particle_position_x,
                         const double *__restrict__ particle_position_y,
                         const double *__restrict__ particle_momentum_z,
                         const short *__restrict__ particle_charge,
                         const double *__restrict__ particle_weight,
                         int *__restrict__ bin_index,
                         int bin_count,
                         const double *__restrict__ invgf_,
                         const int *__restrict__ iold_,
                         const double *__restrict__ deltaold_,
                         double inv_cell_volume,
                         double dx_inv,
                         double dy_inv,
                         double dx_ov_dt,
                         double dy_ov_dt,
                         int    i_domain_begin,
                         int    j_domain_begin,
                         int    nprimy,
                         int    pxr )
{
    naiveCurrentDepositionKernel( Jx,
                                  Jy,
                                  Jz,
                                  Jx_size,
                                  Jy_size,
                                  Jz_size,
                                  particle_position_x,
                                  particle_position_y,
                                  particle_momentum_z,
                                  particle_charge,
                                  particle_weight,
                                  bin_index,
                                  bin_count,
                                  invgf_,
                                  iold_,
                                  deltaold_,
                                  inv_cell_volume,
                                  dx_inv,
                                  dy_inv,
                                  dx_ov_dt,
                                  dy_ov_dt,
                                  i_domain_begin,
                                  j_domain_begin,
                                  nprimy,
                                  pxr );
}

#endif
