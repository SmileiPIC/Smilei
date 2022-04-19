#include "PusherBoris.h"

#include <cmath>
#include <iostream>
#ifdef _GPU
    #include <accelmath.h>
#endif

#include "Particles.h"
#include "Species.h"


PusherBoris::PusherBoris( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherBoris::~PusherBoris()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBoris::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{

    const std::vector<double> *const Epart = &( smpi->dynamics_Epart[ithread] );
    const std::vector<double> *const Bpart = &( smpi->dynamics_Bpart[ithread] );
    double *const __restrict__ invgf       = &( smpi->dynamics_invgf[ithread][0] );

    // double Tx2, Ty2, Tz2;
    // double TxTy, TyTz, TzTx;
    // double alpha;

    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;

    double *const __restrict__ momentum_x = particles.getPtrMomentum( 0 );
    double *const __restrict__ momentum_y = particles.getPtrMomentum( 1 );
    double *const __restrict__ momentum_z = particles.getPtrMomentum( 2 );

    const short *__restrict__ charge = particles.getPtrCharge();

    const int nparts = vecto ? Epart->size() / 3 :
                               particles.last_index.back(); // particles.size()

    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );

#if defined(SMILEI_ACCELERATOR_GPU_OMP)
    const int istart_offset   = istart - ipart_buffer_offset;
    const int particle_number = iend - istart;

    // TODO(Etienne M): Memory ops optimization
    #pragma omp target defaultmap( none ) map( to                                       \
                                               : Ex [istart_offset:particle_number],    \
                                                 Ey [istart_offset:particle_number],    \
                                                 Ez [istart_offset:particle_number],    \
                                                 Bx [istart_offset:particle_number],    \
                                                 By [istart_offset:particle_number],    \
                                                 Bz [istart_offset:particle_number],    \
                                                 invgf [istart_offset:particle_number], \
                                                 charge [istart:particle_number] )      \
        map( tofrom                                                                     \
             : momentum_x [istart:particle_number],                                     \
               momentum_y [istart:particle_number],                                     \
               momentum_z [istart:particle_number],                                     \
               position_x [istart:particle_number],                                     \
               position_y [istart:particle_number],                                     \
               position_z [istart:particle_number] )                                    \
            map( to                                                                     \
                 : istart, iend, ipart_buffer_offset )
    #pragma omp            teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for
#elif defined(_GPU)
    const int istart_offset   = istart - ipart_buffer_offset;
    const int particle_number = iend - istart;

    #pragma acc parallel present(Ex [istart_offset:particle_number],    \
                                 Ey [istart_offset:particle_number],    \
                                 Ez [istart_offset:particle_number],    \
                                 Bx [istart_offset:particle_number],    \
                                 By [istart_offset:particle_number],    \
                                 Bz [istart_offset:particle_number],    \
                                 invgf [istart_offset:particle_number]) \
        deviceptr(position_x,                                           \
                  position_y,                                           \
                  position_z,                                           \
                  momentum_x,                                           \
                  momentum_y,                                           \ 
                  momentum_z,                                           \
                  charge)
    #pragma acc loop gang worker vector
#else
    #pragma omp simd
#endif
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        const int ipart2 = ipart - ipart_buffer_offset;

        const double charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;

        // init Half-acceleration in the electric field
        double pxsm = charge_over_mass_dts2*( Ex[ipart2] );
        double pysm = charge_over_mass_dts2*( Ey[ipart2] );
        double pzsm = charge_over_mass_dts2*( Ez[ipart2] );

        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        const double umx = momentum_x[ipart] + pxsm;
        const double umy = momentum_y[ipart] + pysm;
        const double umz = momentum_z[ipart] + pzsm;

        // Rotation in the magnetic field
        double local_invgf     = charge_over_mass_dts2 / std::sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
        const double Tx        = local_invgf * ( Bx[ipart2] );
        const double Ty        = local_invgf * ( By[ipart2] );
        const double Tz        = local_invgf * ( Bz[ipart2] );
        const double inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );

        pxsm += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* umx  +      2.0*( Tx*Ty+Tz )* umy  +      2.0*( Tz*Tx-Ty )* umz )*inv_det_T;
        pysm += ( 2.0*( Tx*Ty-Tz )* umx  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* umy  +      2.0*( Ty*Tz+Tx )* umz )*inv_det_T;
        pzsm += ( 2.0*( Tz*Tx+Ty )* umx  +      2.0*( Ty*Tz-Tx )* umy  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* umz )*inv_det_T;

        // finalize Half-acceleration in the electric field
        local_invgf = 1. / std::sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        invgf[ipart2] = local_invgf; //1. / std::sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;

        // Move the particle
        local_invgf *= dt;
        // position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        position_x[ipart] += pxsm*local_invgf;
        if( nDim_ > 1 ) {
            position_y[ipart] += pysm*local_invgf;
            if( nDim_ > 2 ) {
                position_z[ipart] += pzsm*local_invgf;
            }
        }
    }

    // if (nDim_>1) {
    //     #pragma omp simd
    //     for( int ipart=istart ; ipart<iend; ipart++ ) {
    //         position_y[ipart] += momentum_y[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
    //     }
    // }
    // 
    // if (nDim_>2) {
    //     #pragma omp simd
    //     for( int ipart=istart ; ipart<iend; ipart++ ) {
    //         position_z[ipart] += momentum_z[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
    //     }
    // }
}
