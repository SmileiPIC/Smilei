#include "PusherBoris.h"

#include <iostream>
#include <cmath>
#ifdef _GPU
#include <accelmath.h>
#endif

#include "Species.h"

#include "Particles.h"

using namespace std;

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

    std::vector<double> * Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> * Bpart = &( smpi->dynamics_Bpart[ithread] );
    double * __restrict__ invgf = &( smpi->dynamics_invgf[ithread][0] );

    double pxsm, pysm, pzsm;
    double local_invgf;
    double umx, umy, umz; //upx, upy, upz;
    double charge_over_mass_dts2;
    double inv_det_T, Tx, Ty, Tz;
    //double Tx2, Ty2, Tz2;
    //double TxTy, TyTz, TzTx;
    //double alpha;

    double * __restrict__ position_x = particles.getPtrPosition(0);
    double * __restrict__ position_y = NULL;
    double * __restrict__ position_z = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
        }
    }
    double * __restrict__ momentum_x = particles.getPtrMomentum(0);
    double * __restrict__ momentum_y = particles.getPtrMomentum(1);
    double * __restrict__ momentum_z = particles.getPtrMomentum(2);

    short * __restrict__ charge = particles.getPtrCharge();

    int nparts;
    if (vecto) {
        nparts = Epart->size()/3;
    } else {
        nparts = particles.size();
    }
    double * Ex = &( ( *Epart )[0*nparts] );
    double * Ey = &( ( *Epart )[1*nparts] );
    double * Ez = &( ( *Epart )[2*nparts] );
    double * Bx = &( ( *Bpart )[0*nparts] );
    double * By = &( ( *Bpart )[1*nparts] );
    double * Bz = &( ( *Bpart )[2*nparts] );

    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        int ipart2 = ipart - ipart_buffer_offset;

        charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;

        // init Half-acceleration in the electric field
        pxsm = charge_over_mass_dts2*( *( Ex+ipart2 ) );
        pysm = charge_over_mass_dts2*( *( Ey+ipart2 ) );
        pzsm = charge_over_mass_dts2*( *( Ez+ipart2) );

        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        umx = momentum_x[ipart] + pxsm;
        umy = momentum_y[ipart] + pysm;
        umz = momentum_z[ipart] + pzsm;

        // Rotation in the magnetic field
        local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
        Tx    = local_invgf * ( *( Bx+ipart2 ) );
        Ty    = local_invgf * ( *( By+ipart2 ) );
        Tz    = local_invgf * ( *( Bz+ipart2 ) );
        inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );

        pxsm += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* umx  +      2.0*( Tx*Ty+Tz )* umy  +      2.0*( Tz*Tx-Ty )* umz )*inv_det_T;
        pysm += ( 2.0*( Tx*Ty-Tz )* umx  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* umy  +      2.0*( Ty*Tz+Tx )* umz )*inv_det_T;
        pzsm += ( 2.0*( Tz*Tx+Ty )* umx  +      2.0*( Ty*Tz-Tx )* umy  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* umz )*inv_det_T;

        // finalize Half-acceleration in the electric field
        //local_invgf = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        invgf[ipart2] = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;

        // Move the particle
        //local_invgf *= dt;
        position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart2];
        // if (nDim_>1) {
        //     position_y[ipart] += pysm*local_invgf;
        //     if (nDim_>2) {
        //         position_z[ipart] += pzsm*local_invgf;
        //     }
        // }
    }

    if (nDim_>1) {
        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            position_y[ipart] += momentum_y[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
        }
    }

    if (nDim_>2) {
        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            position_z[ipart] += momentum_z[ipart]*invgf[ipart-ipart_buffer_offset]*dt;
        }
    }
}
