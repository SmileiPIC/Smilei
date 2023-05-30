#include "PusherBorisBTIS3.h"

#include <iostream>
#include <cmath>
#ifdef _GPU
#include <accelmath.h>
#endif

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherBorisBTIS3::PusherBorisBTIS3( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherBorisBTIS3::~PusherBorisBTIS3()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBorisBTIS3::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{

    const int nparts = vecto ? smpi->dynamics_Epart[ithread].size() / 3 :
                               particles.last_index.back(); // particles.size()

    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;
    
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);

    const short *const __restrict__ charge = particles.getPtrCharge();

    double *const __restrict__ invgf = &( smpi->dynamics_invgf[ithread][0] );
    const double *const __restrict__ Ex      = &( ( smpi->dynamics_Epart[ithread] )[0*nparts] );
    const double *const __restrict__ Ey      = &( ( smpi->dynamics_Epart[ithread] )[1*nparts] );
    const double *const __restrict__ Ez      = &( ( smpi->dynamics_Epart[ithread] )[2*nparts] );
    const double *const __restrict__ Bx      = &( ( smpi->dynamics_Bpart[ithread] )[0*nparts] );
    const double *const __restrict__ By      = &( ( smpi->dynamics_Bpart[ithread] )[1*nparts] );
    const double *const __restrict__ Bz      = &( ( smpi->dynamics_Bpart[ithread] )[2*nparts] );
    const double *const __restrict__ ByBTIS3 = &( ( smpi->dynamics_Bpart_yBTIS3[ithread] )[0*nparts] );
    const double *const __restrict__ BzBTIS3 = &( ( smpi->dynamics_Bpart_zBTIS3[ithread] )[0*nparts] );

    #pragma omp simd
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
        const double Tx        = local_invgf * ( Bx     [ipart2] );
        const double Ty        = local_invgf * ( By     [ipart2] );
        const double Tz        = local_invgf * ( Bz     [ipart2] );
        const double TyBTIS3   = local_invgf * ( ByBTIS3[ipart2] );
        const double TzBTIS3   = local_invgf * ( BzBTIS3[ipart2] );
        // const double inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz ); // classic Boris
        const double inv_det_T = 1./( 1.+Tx*Tx+Ty*TyBTIS3+Tz*TzBTIS3+TzBTIS3*Tx*Ty-Tx*TyBTIS3*Tz ); // Boris BTIS3
 
        // classic Boris

        // pxsm += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* umx  +      2.0*( Tx*Ty+Tz )* umy  +      2.0*( Tz*Tx-Ty )* umz )*inv_det_T;
        // pysm += ( 2.0*( Tx*Ty-Tz )* umx  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* umy  +      2.0*( Ty*Tz+Tx )* umz )*inv_det_T;
        // pzsm += ( 2.0*( Tz*Tx+Ty )* umx  +      2.0*( Ty*Tz-Tx )* umy  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* umz )*inv_det_T;

        // or also, for the sake of comparison with BTIS3 modification:

        // pxsm += ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )  * umx * inv_det_T; 
        // pxsm += 2.0 * ( Tx*Ty+Tz )         * umy * inv_det_T;  
        // pxsm += 2.0 * ( Tz*Tx-Ty )         * umz * inv_det_T;
        // 
        // pysm += 2.0 * ( Tx*Ty-Tz )         * umx * inv_det_T;
        // pysm += ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )  * umy * inv_det_T;
        // pysm += 2.0 * ( Ty*Tz+Tx )         * umz * inv_det_T;
        // 
        // pzsm += 2.0 * ( Tz*Tx+Ty )         * umx * inv_det_T;
        // pzsm += 2.0 * ( Ty*Tz-Tx )         * umy * inv_det_T;
        // pzsm += ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )  * umz * inv_det_T;

        // borisBTIS3
        pxsm += ( 1.0 + Tx*Tx      -TzBTIS3*(Tx*Ty + Tz)  + TyBTIS3*(Tx*Tz - Ty)  )     * umx * inv_det_T; 
        pxsm += 2.0 * (Tx*Ty + Tz)                                                      * umy * inv_det_T;  
        pxsm += 2.0 * (Tz*Tx - Ty)                                                      * umz * inv_det_T;

        pysm += 2.0 * ( Tx*TyBTIS3-TzBTIS3 )                                            * umx * inv_det_T;
        pysm += ( 1.0 + Ty*TyBTIS3 - Tx*(TzBTIS3*Ty + Tx) - Tz*(TzBTIS3 - Tx*TyBTIS3))  * umy * inv_det_T;
        pysm += 2.0 * (TzBTIS3*Ty + Tx)                                                 * umz * inv_det_T;

        pzsm += 2.0 * (TzBTIS3*Tx + TyBTIS3)                                            * umx * inv_det_T;
        pzsm += 2.0 * (-Tx + TyBTIS3*Tz)                                                * umy * inv_det_T;
        pzsm += ( 1.0+ TzBTIS3*Tz - Tx*(Tx - TyBTIS3*Tz) - Ty*(TzBTIS3*Tx + TyBTIS3)  ) * umz * inv_det_T;
        

        // finalize Half-acceleration in the electric field
        local_invgf = 1. / std::sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        invgf[ipart2] = local_invgf; //1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;

        // Move the particle
        local_invgf *= dt;
        //position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        position_x[ipart] += pxsm*local_invgf;
        if (nDim_>1) {
            position_y[ipart] += pysm*local_invgf;
            if (nDim_>2) {
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
