/*! @file PusherVay.cpp

 @brief PusherVay.cpp  generic class for the particle pusher of J.L. Vay.

 #details The description of the J.L. Vay pusher can be found in this reference:
          http://dx.doi.org/10.1063/1.2837054

 @date 2017-04-10
*/
#include "PusherVay.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherVay::PusherVay( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherVay::~PusherVay()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Vay) scheme
***********************************************************************/

void PusherVay::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    double charge_over_mass_dts2;
    double upx, upy, upz, us2;
    double alpha, s, T2 ;
    double Tx, Ty, Tz;
    double pxsm, pysm, pzsm;
    // Only useful for the second method
    //double Tx2, Ty2, Tz2;
    //double TxTy, TyTz, TzTx;
    
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
        
        // ____________________________________________
        // Part I: Computation of uprime
        
        // For unknown reason, this has to be computed again
        ( *invgf )[ipart] = 1./sqrt( 1.0 + momentum[0][ipart]*momentum[0][ipart]
                                     + momentum[1][ipart]*momentum[1][ipart]
                                     + momentum[2][ipart]*momentum[2][ipart] );
                                     
        // Add Electric field
        upx = momentum[0][ipart] + 2.*charge_over_mass_dts2*( *( Ex+ipart ) );
        upy = momentum[1][ipart] + 2.*charge_over_mass_dts2*( *( Ey+ipart ) );
        upz = momentum[2][ipart] + 2.*charge_over_mass_dts2*( *( Ez+ipart ) );
        
        // Add magnetic field
        Tx  = charge_over_mass_dts2* ( *( Bx+ipart ) );
        Ty  = charge_over_mass_dts2* ( *( By+ipart ) );
        Tz  = charge_over_mass_dts2* ( *( Bz+ipart ) );
        
        upx += ( *invgf )[ipart]*( momentum[1][ipart]*Tz - momentum[2][ipart]*Ty );
        upy += ( *invgf )[ipart]*( momentum[2][ipart]*Tx - momentum[0][ipart]*Tz );
        upz += ( *invgf )[ipart]*( momentum[0][ipart]*Ty - momentum[1][ipart]*Tx );
        
        // alpha is gamma^2
        alpha = 1.0 + upx*upx + upy*upy + upz*upz;
        T2    = Tx*Tx + Ty*Ty + Tz*Tz;
        
        // ___________________________________________
        // Part II: Computation of Gamma^{i+1}
        
        // s is sigma
        s     = alpha - T2;
        us2   = pow( upx*Tx + upy*Ty + upz*Tz, 2.0 );
        
        // alpha becomes 1/gamma^{i+1}
        alpha = 1.0/sqrt( 0.5*( s + sqrt( s*s + 4.0*( T2 + us2 ) ) ) );
        
        Tx *= alpha;
        Ty *= alpha;
        Tz *= alpha;
        
        s = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );
        alpha   = upx*Tx + upy*Ty + upz*Tz;
        
        pxsm = s*( upx + alpha*Tx + Tz*upy - Ty*upz );
        pysm = s*( upy + alpha*Ty + Tx*upz - Tz*upx );
        pzsm = s*( upz + alpha*Tz + Ty*upx - Tx*upy );
        
        // Second way of doing it like in the Boris pusher
        //Tx2   = Tx*Tx;
        //Ty2   = Ty*Ty;
        //Tz2   = Tz*Tz;
        
        //TxTy  = Tx*Ty;
        //TyTz  = Ty*Tz;
        //TzTx  = Tz*Tx;
        
        //pxsm = ((1.0+Tx2)* upx  + (TxTy+Tz)* upy + (TzTx-Ty)* upz)*s;
        //pysm = ((TxTy-Tz)* upx  + (1.0+Ty2)* upy + (TyTz+Tx)* upz)*s;
        //pzsm = ((TzTx+Ty)* upx  + (TyTz-Tx)* upy + (1.0+Tz2)* upz)*s;
        
        // Inverse Gamma factor
        ( *invgf )[ipart] = 1.0 / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        
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
        
        particles.cell_keys.resize( nparts );
        cell_keys = &( particles.cell_keys[0] );
        #pragma omp simd
        for( int ipart=0 ; ipart<nparts; ipart++ ) {
        
            for( int i = 0 ; i<nDim_ ; i++ ) {
                cell_keys[ipart] *= nspace[i];
                cell_keys[ipart] += round( ( position[i][ipart]-min_loc_vec[i] ) * dx_inv_[i] );
            }
        }
    }
    
    
    
}
