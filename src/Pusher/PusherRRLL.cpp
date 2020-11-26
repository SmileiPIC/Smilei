#include "PusherRRLL.h"

#include <iostream>
#include <cmath>

#include "Species.h"
#include "Particles.h"

using namespace std;

PusherRRLL::PusherRRLL( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherRRLL::~PusherRRLL()
{
}

/****************************************************************************
    Lorentz Force -- leap-frog (Boris) scheme + classical rad. reaction force
*****************************************************************************/
void PusherRRLL::operator()( Particles &particles, int ipart, LocalFields Epart, LocalFields Bpart, double &invgf )
{
    // Declaration of local variables
    // ------------------------------
    
    double charge_over_mass_ = static_cast<double>( particles.charge( ipart ) )*one_over_mass_;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    
    //DEBUG(5, "\tPush particle"<< particles.position(0, ipart) );
    
    // --------------------------------------
    // SOLVE THE PARTICLE EQUATION OF MOTIONS
    // --------------------------------------
    
    // Half-acceleration in the electric field
    umx = particles.momentum( 0, ipart ) + charge_over_mass_*Epart.x*dts2;
    umy = particles.momentum( 1, ipart ) + charge_over_mass_*Epart.y*dts2;
    umz = particles.momentum( 2, ipart ) + charge_over_mass_*Epart.z*dts2;
    invgf  = 1. / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
    
    // Rotation in the magnetic field
    alpha = charge_over_mass_*dts2*invgf;
    Tx    = alpha * Bpart.x;
    Ty    = alpha * Bpart.y;
    Tz    = alpha * Bpart.z;
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
    
    // Half-acceleration in the electric field
    pxsm = upx + charge_over_mass_*Epart.x*dts2;
    pysm = upy + charge_over_mass_*Epart.y*dts2;
    pzsm = upz + charge_over_mass_*Epart.z*dts2;
    invgf = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
    
    particles.momentum( 0, ipart ) = pxsm;
    particles.momentum( 1, ipart ) = pysm;
    particles.momentum( 2, ipart ) = pzsm;
    
    // Move the particle
    for( int i = 0 ; i<nDim_ ; i++ ) {
        particles.position( i, ipart )     += dt*particles.momentum( i, ipart )*invgf;
    }
    
    // COMPUTE Chi
    particles.chi( ipart )=0.5;
    //DEBUG(5, "\t END "<< particles.position(0, ipart) );
    
}
void PusherRRLL::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref )
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
    double charge_over_mass_ ;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double local_invgf;
    
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        //(*this)(particles, iPart, (*Epart)[iPart], (*Bpart)[iPart] , (*invgf)[iPart]);
        charge_over_mass_ = static_cast<double>( particles.charge( ipart ) )*one_over_mass_;
        // Half-acceleration in the electric field
        umx = particles.momentum( 0, ipart ) + charge_over_mass_*( *( Ex+ipart ) )*dts2;
        umy = particles.momentum( 1, ipart ) + charge_over_mass_*( *( Ey+ipart ) )*dts2;
        umz = particles.momentum( 2, ipart ) + charge_over_mass_*( *( Ez+ipart ) )*dts2;
        local_invgf  = 1. / sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
        
        // Rotation in the magnetic field
        alpha = charge_over_mass_*dts2*local_invgf;
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
        
        // Half-acceleration in the electric field
        pxsm = upx + charge_over_mass_*( *( Ex+ipart ) )*dts2;
        pysm = upy + charge_over_mass_*( *( Ey+ipart ) )*dts2;
        pzsm = upz + charge_over_mass_*( *( Ez+ipart ) )*dts2;
        ( *invgf )[ipart] = 1. / sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
        
        particles.momentum( 0, ipart ) = pxsm;
        particles.momentum( 1, ipart ) = pysm;
        particles.momentum( 2, ipart ) = pzsm;
        
        // Move the particle
        for( int i = 0 ; i<nDim_ ; i++ ) {
            particles.position( i, ipart )     += dt*particles.momentum( i, ipart )*( *invgf )[ipart];
        }
        
        // COMPUTE Chi
        particles.chi( ipart )=0.5;
        //DEBUG(5, "\t END "<< particles.position(0, ipart) );
    }
    
    if( vecto ) {
        double *position[3];
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i] =  &( particles.position( i, 0 ) );
        }
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
