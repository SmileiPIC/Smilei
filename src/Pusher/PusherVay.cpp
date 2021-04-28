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

void PusherVay::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    double *invgf = &( smpi->dynamics_invgf[ithread][0] );
    
    double charge_over_mass_dts2;
    double upx, upy, upz, us2;
    double alpha, s, T2 ;
    double Tx, Ty, Tz;
    double pxsm, pysm, pzsm;
    
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

    short *charge = particles.getPtrCharge();
    
    int nparts;
    if (vecto) {
        nparts = Epart->size()/3;
    } else {
        //nparts = particles.size();
        nparts = particles.last_index.back();
    }
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    #ifndef _GPU
        #pragma omp simd
    #else
        int np = iend-istart;
        #pragma acc parallel present(Ex[istart:np],Ey[istart:np],Ez[istart:np],Bx[istart:np],By[istart:np],Bz[istart:np],invgf[0:nparts]) deviceptr(position_x,position_y,position_z,momentum_x,momentum_y,momentum_z,charge)
        #pragma acc loop gang worker vector
    #endif
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass_dts2 = ( double )( charge[ipart] )*one_over_mass_*dts2;
        
        // ____________________________________________
        // Part I: Computation of uprime
        
        // For unknown reason, this has to be computed again
        invgf[ipart-ipart_buffer_offset] = 1./sqrt( 1.0 + momentum_x[ipart]*momentum_x[ipart]
                                     + momentum_y[ipart]*momentum_y[ipart]
                                     + momentum_z[ipart]*momentum_z[ipart] );
                                     
        // Add Electric field
        upx = momentum_x[ipart] + 2.*charge_over_mass_dts2*( *( Ex+ipart-ipart_buffer_offset ) );
        upy = momentum_y[ipart] + 2.*charge_over_mass_dts2*( *( Ey+ipart-ipart_buffer_offset ) );
        upz = momentum_z[ipart] + 2.*charge_over_mass_dts2*( *( Ez+ipart-ipart_buffer_offset ) );
        
        // Add magnetic field
        Tx  = charge_over_mass_dts2* ( *( Bx+ipart-ipart_buffer_offset ) );
        Ty  = charge_over_mass_dts2* ( *( By+ipart-ipart_buffer_offset ) );
        Tz  = charge_over_mass_dts2* ( *( Bz+ipart-ipart_buffer_offset ) );
        
        upx += invgf [ipart-ipart_buffer_offset]*( momentum_y[ipart]*Tz - momentum_z[ipart]*Ty );
        upy += invgf [ipart-ipart_buffer_offset]*( momentum_z[ipart]*Tx - momentum_x[ipart]*Tz );
        upz += invgf [ipart-ipart_buffer_offset]*( momentum_x[ipart]*Ty - momentum_y[ipart]*Tx );
        
        // alpha is gamma^2
        alpha = 1.0 + upx*upx + upy*upy + upz*upz;
        T2    = Tx*Tx + Ty*Ty + Tz*Tz;
        
        // ___________________________________________
        // Part II: Computation of Gamma^{i+1}
        
        // s is sigma
        s     = alpha - T2;
        us2   = upx*Tx + upy*Ty + upz*Tz;
        us2   = us2*us2;
        
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
        invgf[ipart-ipart_buffer_offset] = 1.0 / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        
        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;
        
        // Move the particle
        position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart-ipart_buffer_offset];
        if (nDim_>1) {
            position_y[ipart] += dt*momentum_y[ipart]*invgf[ipart-ipart_buffer_offset];
            if (nDim_>2) {
                position_z[ipart] += dt*momentum_z[ipart]*invgf[ipart-ipart_buffer_offset];
            }
        }
        
    }
    
    // if( vecto ) {
    //     int *cell_keys;
    //
    //     particles.cell_keys.resize( nparts );
    //     cell_keys = &( particles.cell_keys[0] );
    //     #pragma omp simd
    //     for( int ipart=0 ; ipart<nparts; ipart++ ) {
    //
    //         for( int i = 0 ; i<nDim_ ; i++ ) {
    //             cell_keys[ipart] *= nspace[i];
    //             cell_keys[ipart] += round( ( position[i][ipart]-min_loc_vec[i] ) * dx_inv_[i] );
    //         }
    //     }
    // }
    
}
