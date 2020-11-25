#include "PusherBoris.h"

#include <iostream>
#include <cmath>

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

void PusherBoris::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    double charge_over_mass_dts2;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double local_invgf;
    
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, 0 ) );
    }
    double *position[3];
    for( int i = 0 ; i<nDim_ ; i++ ) {
        position[i] =  &( particles.position( i, 0 ) );
    }
    
    short *charge = &( particles.charge( 0 ) );
    
    int nparts = particles.size();
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;
        
        // init Half-acceleration in the electric field
        pxsm = charge_over_mass_dts2*( *( Ex+ipart ) );
        pysm = charge_over_mass_dts2*( *( Ey+ipart ) );
        pzsm = charge_over_mass_dts2*( *( Ez+ipart ) );
        
        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        umx = momentum[0][ipart] + pxsm;
        umy = momentum[1][ipart] + pysm;
        umz = momentum[2][ipart] + pzsm;
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
        ( *invgf )[ipart] = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        
        momentum[0][ipart] = pxsm;
        momentum[1][ipart] = pysm;
        momentum[2][ipart] = pzsm;
        
        // Move the particle
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i][ipart]     += dt*momentum[i][ipart]*( *invgf )[ipart];
        }
        
    }
}
