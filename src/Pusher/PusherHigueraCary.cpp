/*! @file PusherHigueraCary.cpp

  @brief PusherHigueraCary.cpp  generic class for the particle pusher of A.V. Higuera and J.R. Cari.

  @details See article https://arxiv.org/abs/1701.05605

  @date 2017-03-31
*/

#include "PusherHigueraCary.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherHigueraCary::PusherHigueraCary( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherHigueraCary::~PusherHigueraCary()
{
}

/***********************************************************************
  Lorentz Force -- leap-frog (HigueraCary) scheme
 ***********************************************************************/

void PusherHigueraCary::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    
    int nparts = particles.size();
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    double charge_over_mass_dts2;
    double umx, umy, umz, upx, upy, upz, gfm2;
    double beta2, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
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
        
        // Intermediate gamma factor: only this part differs from the Boris scheme
        // Square Gamma factor from um
        gfm2 = ( 1.0 + umx*umx + umy*umy + umz*umz );
        
        // Equivalent of betax,betay,betaz in the paper
        Tx    = charge_over_mass_dts2 * ( *( Bx+ipart ) );
        Ty    = charge_over_mass_dts2 * ( *( By+ipart ) );
        Tz    = charge_over_mass_dts2 * ( *( Bz+ipart ) );
        
        // beta**2
        beta2 = Tx*Tx + Ty*Ty + Tz*Tz;
        
        // Equivalent of 1/\gamma_{new} in the paper
        local_invgf = 1./sqrt( 0.5*( gfm2 - beta2 +
                                     sqrt( pow( gfm2 - beta2, 2 ) + 4.0*( beta2 + pow( Tx*umx + Ty*umy + Tz*umz, 2 ) ) ) ) );
                                     
        // Rotation in the magnetic field
        Tx    *= local_invgf;
        Ty    *= local_invgf;
        Tz    *= local_invgf;
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
        
        // final gamma factor
        ( *invgf )[ipart] = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        
        momentum[0][ipart] = pxsm;
        momentum[1][ipart] = pysm;
        momentum[2][ipart] = pzsm;
        
        // Move the particle
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i][ipart]     += dt*momentum[i][ipart]*( *invgf )[ipart];
        }
        
    }
    
    if( vecto ) {
        int *cell_keys;
        particles.cell_keys.resize( iend-istart );
        cell_keys = &( particles.cell_keys[0] );
        
        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
        
            for( int i = 0 ; i<nDim_ ; i++ ) {
                cell_keys[ipart] *= nspace[i];
                cell_keys[ipart] += round( ( position[i][ipart]-min_loc_vec[i] ) * dx_inv_[i] );
            }
            
        }
    }
}
