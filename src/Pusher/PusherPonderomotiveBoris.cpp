#include "PusherPonderomotiveBoris.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;
// Pushes only momentum of particles interacting with envelope, not their position
PusherPonderomotiveBoris::PusherPonderomotiveBoris( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherPonderomotiveBoris::~PusherPonderomotiveBoris()
{
}

/**************************************************************************
    Lorentz Force + Ponderomotive force -- leap-frog (Boris-style) scheme, momentum advance
**************************************************************************/

void PusherPonderomotiveBoris::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{
    std::vector<double> *Epart       = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart       = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *GradPhipart = &( smpi->dynamics_GradPHIpart[ithread] );
    double *dynamics_inv_gamma_ponderomotive = &( smpi->dynamics_inv_gamma_ponderomotive[ithread][0] );
    
    double charge_over_mass_dts2, charge_sq_over_mass_sq_dts4;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    //double one_ov_gamma_ponderomotive;
    
    double *const __restrict__ momentum_x = particles.getPtrMomentum(0);
    double *const __restrict__ momentum_y = particles.getPtrMomentum(1);
    double *const __restrict__ momentum_z = particles.getPtrMomentum(2);
    
    short *charge = particles.getPtrCharge();
    
    const int nparts = vecto ? Epart->size() / 3 :
                               particles.size(); // particles.size()
    
    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );
    const double *const __restrict__ GradPhix = &( ( *GradPhipart )[0*nparts] );
    const double *const __restrict__ GradPhiy = &( ( *GradPhipart )[1*nparts] );
    const double *const __restrict__ GradPhiz = &( ( *GradPhipart )[2*nparts] );
    //double *inv_gamma_ponderomotive = &( ( *dynamics_inv_gamma_ponderomotive )[0*nparts] );
    
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
    
        charge_over_mass_dts2    = ( double )( charge[ipart] )*one_over_mass_*dts2;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_sq_dts4 = ( double )( charge[ipart] )*( double )( charge[ipart] )*one_over_mass_*one_over_mass_*dts4;
        
        // ponderomotive gamma buffered from susceptibility
        // one_ov_gamma_ponderomotive = dynamics_inv_gamma_ponderomotive[ipart-ipart_buffer_offset];
        
        // init Half-acceleration in the electric field and ponderomotive force
        pxsm = charge_over_mass_dts2 * ( Ex[ipart-ipart_buffer_offset] ) - charge_sq_over_mass_sq_dts4 * ( GradPhix[ipart-ipart_buffer_offset] ) * dynamics_inv_gamma_ponderomotive[ipart-ipart_buffer_offset] ;
        pysm = charge_over_mass_dts2 * ( Ey[ipart-ipart_buffer_offset] ) - charge_sq_over_mass_sq_dts4 * ( GradPhiy[ipart-ipart_buffer_offset] ) * dynamics_inv_gamma_ponderomotive[ipart-ipart_buffer_offset] ;
        pzsm = charge_over_mass_dts2 * ( Ez[ipart-ipart_buffer_offset] ) - charge_sq_over_mass_sq_dts4 * ( GradPhiz[ipart-ipart_buffer_offset] ) * dynamics_inv_gamma_ponderomotive[ipart-ipart_buffer_offset] ;
        
        umx = momentum_x[ipart] + pxsm;
        umy = momentum_y[ipart] + pysm;
        umz = momentum_z[ipart] + pzsm;
        
        // Rotation in the magnetic field, using updated gamma ponderomotive
        alpha = charge_over_mass_dts2 * dynamics_inv_gamma_ponderomotive[ipart-ipart_buffer_offset];
        Tx    = alpha * ( Bx[ipart-ipart_buffer_offset] );
        Ty    = alpha * ( By[ipart-ipart_buffer_offset] );
        Tz    = alpha * ( Bz[ipart-ipart_buffer_offset] );
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
        
        // finalize Half-acceleration in the electric field and ponderomotive force
        pxsm += upx;
        pysm += upy;
        pzsm += upz;
        
        momentum_x[ipart] = pxsm;
        momentum_y[ipart] = pysm;
        momentum_z[ipart] = pzsm;
        
    }
}
