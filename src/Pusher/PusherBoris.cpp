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
    
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    double *invgf = &( smpi->dynamics_invgf[ithread][0] );

    double* position_x = particles.getPtrPosition(0);
    double* position_y = NULL;
    double* position_z = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
        }
    }
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    
    short *charge = &( particles.charge( 0 ) );
    
    int nparts;
    if (vecto) {
        nparts = Epart->size()/3;
    } else {
        nparts = particles.size();
    }
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );

    if (vecto) {

        double charge_over_mass_dts2;
        double inv_det_T, Tx, Ty, Tz;
        double local_invgf;
        //int IX;

        //int* cell_keys;

        //particles.cell_keys.resize(nparts);
        //cell_keys = &( particles.cell_keys[0]);

        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            double psm[3], um[3];
            
            charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;
            
            // init Half-acceleration in the electric field
            psm[0] = charge_over_mass_dts2*( *( Ex+ipart-ipart_buffer_offset ) );
            psm[1] = charge_over_mass_dts2*( *( Ey+ipart-ipart_buffer_offset ) );
            psm[2] = charge_over_mass_dts2*( *( Ez+ipart-ipart_buffer_offset ) );
            
            //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
            um[0] = momentum_x[ipart] + psm[0];
            um[1] = momentum_y[ipart] + psm[1];
            um[2] = momentum_z[ipart] + psm[2];
        
            // Rotation in the magnetic field
            local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + um[0]*um[0] + um[1]*um[1] + um[2]*um[2] );
            Tx    = local_invgf * ( *( Bx+ipart-ipart_buffer_offset ) );
            Ty    = local_invgf * ( *( By+ipart-ipart_buffer_offset ) );
            Tz    = local_invgf * ( *( Bz+ipart-ipart_buffer_offset ) );
            inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );
            
            psm[0] += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* um[0]  +      2.0*( Tx*Ty+Tz )* um[1]  +      2.0*( Tz*Tx-Ty )* um[2] )*inv_det_T;
            psm[1] += ( 2.0*( Tx*Ty-Tz )* um[0]  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* um[1]  +      2.0*( Ty*Tz+Tx )* um[2] )*inv_det_T;
            psm[2] += ( 2.0*( Tz*Tx+Ty )* um[0]  +      2.0*( Ty*Tz-Tx )* um[1]  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* um[2] )*inv_det_T;
            
            // finalize Half-acceleration in the electric field
            local_invgf = 1. / sqrt( 1.0 + psm[0]*psm[0] + psm[1]*psm[1] + psm[2]*psm[2] );
            invgf[ipart-ipart_buffer_offset] = local_invgf;
            
            momentum_x[ipart] = psm[0];
            momentum_y[ipart] = psm[1];
            momentum_z[ipart] = psm[2];
            
            // Move the particle
            local_invgf *= dt;
            position_x[ipart] += psm[0]*local_invgf;
            if (nDim_>1) {
                position_y[ipart] += psm[1]*local_invgf;
                if (nDim_>2) {
                    position_z[ipart] += psm[2]*local_invgf;
                }
            }
        }

        
    } else {
        
        double charge_over_mass_dts2;
        double umx, umy, umz, upx, upy, upz;
        double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
        double TxTy, TyTz, TzTx;
        double pxsm, pysm, pzsm;
        double local_invgf;
    
        
        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
            charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;
            
            // init Half-acceleration in the electric field
            pxsm = charge_over_mass_dts2*( *( Ex+ipart ) );
            pysm = charge_over_mass_dts2*( *( Ey+ipart ) );
            pzsm = charge_over_mass_dts2*( *( Ez+ipart ) );
            
            //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
            umx = momentum_x[ipart] + pxsm;
            umy = momentum_y[ipart] + pysm;
            umz = momentum_z[ipart] + pzsm;
            local_invgf = 1. / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
            
            // Rotation in the magnetic field
            alpha = charge_over_mass_dts2*local_invgf;
            Tx    = alpha * ( *( Bx+ipart ) );
            Ty    = alpha * ( *( By+ipart ) );
            Tz    = alpha * ( *( Bz+ipart ) );
            Tx2   = Tx*Tx;
            Ty2   = Ty*Ty;
            Tz2   = Tz*Tz;
            TxTy  = Tx*Ty;
            TyTz  = Ty*Tz;
            TzTx  = Tz*Tx;
            inv_det_T = 1.0/( 1.0+Tx2+Ty2+Tz2 );
            
            upx = ( ( 1.0+Tx2-Ty2-Tz2 )* umx  +      2.0*( TxTy+Tz )* umy  +      2.0*( TzTx-Ty )* umz )*inv_det_T;
            upy = ( 2.0*( TxTy-Tz )* umx  + ( 1.0-Tx2+Ty2-Tz2 )* umy  +      2.0*( TyTz+Tx )* umz )*inv_det_T;
            upz = ( 2.0*( TzTx+Ty )* umx  +      2.0*( TyTz-Tx )* umy  + ( 1.0-Tx2-Ty2+Tz2 )* umz )*inv_det_T;
            
            // finalize Half-acceleration in the electric field
            pxsm += upx;
            pysm += upy;
            pzsm += upz;
            invgf[ipart] = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
            
            momentum_x[ipart] = pxsm;
            momentum_y[ipart] = pysm;
            momentum_z[ipart] = pzsm;
            
            // Move the particle
            position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart];
            if (nDim_>1) {
                position_y[ipart] += dt*momentum_y[ipart]*invgf[ipart];
                if (nDim_>2) {
                    position_z[ipart] += dt*momentum_z[ipart]*invgf[ipart];
                }
            }
        }
    }
}
