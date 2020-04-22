#include "InterpolatorAM1Order.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include "ElectroMagn.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include <complex>
#include "dcomplex.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for InterpolatorAM1Order
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorAM1Order::InterpolatorAM1Order( Params &params, Patch *patch ) : InterpolatorAM( params, patch )
{

    dl_inv_ = 1.0/params.cell_length[0];
    dr_inv_ = 1.0/params.cell_length[1];
    nmodes = params.nmodes;
    dr =  params.cell_length[1];
}

// ---------------------------------------------------------------------------------------------------------------------
// 1st Order Interpolation of the fields at a the particle position (2 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void InterpolatorAM1Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{

    //Treat mode 0 first
    
    // Static cast of the electromagnetic fields
    cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
    cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
    cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
    cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
    cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
    cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
    
    
    // Normalized particle position
    // Spectral grids are shifted by a half cell length along r.
    double xpn = particles.position( 0, ipart ) * dl_inv_ ;
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * dr_inv_ - 0.5 ; //-0.5 because of cells being shifted by dr/2
    exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ; //exp(-i theta)
    complex<double> exp_mm_theta = 1. ;                                                          //exp(-i m theta)
    // Calculate coeffs
    coeffs( xpn, rpn );
 
    // Interpolation of El^(p,p)
    *( ELoc+0*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[0], El, ip_, jp_ ) );
    // Interpolation of Er^(p,p) zero on axis
    *( ELoc+1*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Er, ip_, jp_ ) );
    // Interpolation of Et^(p,p) zero on axis
    *( ELoc+2*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Et, ip_, jp_ ) );
    // Interpolation of Bl^(p,p)
    *( BLoc+0*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[0], Bl, ip_, jp_ ) );
    // Interpolation of Br^(p,p) zero on axis
    *( BLoc+1*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Br, ip_, jp_ ) );
    // Interpolation of Bt^(p,p) zero on axis
    *( BLoc+2*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Bt, ip_, jp_ ) );
    
    for( unsigned int imode = 1; imode < nmodes ; imode++ ) {
        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
        
        exp_mm_theta *= exp_m_theta ;
        //int zero_on_axis_l_rho  =  (imode != 0)*2;//=2 if El, Bl, Jl or rho is zero on axis and 0 if they are constant.
        int zero_on_axis_rt = (imode != 1)*2;     //=2 if Er, Et, Br, Bt, Jr, Jt is zero on axis and 0 if they are constant.

        *( ELoc+0*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[2], El, ip_, jp_ )* exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Er, ip_, jp_ )* exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Et, ip_, jp_ )* exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[2], Bl, ip_, jp_ )* exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Br, ip_, jp_ )* exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Bt, ip_, jp_ )* exp_mm_theta ) ;
    }
    
    //Translate field into the cartesian y,z coordinates
    double delta2 = std::real( exp_m_theta ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta ) * *( ELoc+1*nparts ) + std::real( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta ) * *( BLoc+1*nparts ) + std::real( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;
    
} // END InterpolatorAM1Order

//Interpolator for probes
void InterpolatorAM1Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    int ipart = *istart;
    
    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );
    
    // Interpolate E, B
    cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
    cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
    cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
    cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
    cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
    cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
    cField2D *Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[0];
    cField2D *Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[0];
    cField2D *Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[0];
    cField2D *Rho= ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[0];
    
    // Normalized particle position
    double xpn = particles.position( 0, ipart ) * dl_inv_;
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * dr_inv_ - 0.5 ; //-0.5 because of cells being shifted by dr/2
    complex<double> exp_mm_theta = 1. ;
    
    // Calculate coeffs
    coeffs( xpn, rpn );
 
    int nparts( particles.size() );
    
    // Interpolation of El^(p,p)
    *( ELoc+0*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[0], El, ip_, jp_ ) );
    // Interpolation of Er^(p,p) zero on axis
    *( ELoc+1*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Er, ip_, jp_ ) );
    // Interpolation of Et^(p,p) zero on axis
    *( ELoc+2*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Et, ip_, jp_ ) );
    // Interpolation of Bl^(p,p)
    *( BLoc+0*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[0], Bl, ip_, jp_ ) );
    // Interpolation of Br^(p,p) zero on axis
    *( BLoc+1*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Br, ip_, jp_ ) );
    // Interpolation of Bt^(p,p) zero on axis
    *( BLoc+2*nparts ) = std::real( compute( &coeffxp_[0], &coeffyp_[2], Bt, ip_, jp_ ) );
    // Interpolation of Jl^(p,p,p)
    JLoc->x = std::real( compute( &coeffxp_[0], &coeffyp_[0], Jl, ip_, jp_ ) );
    // Interpolation of Jr^(p,p,p) zero on axis
    JLoc->y = std::real( compute( &coeffxp_[0], &coeffyp_[2], Jr, ip_, jp_ ) );
    // Interpolation of Jt^(p,p,p) zero on axis
    JLoc->z = std::real( compute( &coeffxp_[0], &coeffyp_[2], Jt, ip_, jp_ ) );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = std::real( compute( &coeffxp_[0], &coeffyp_[0], Rho, ip_, jp_ ) );
   
    if (r > 0){ 
        exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
    } else {
        exp_m_theta = 1. ;
    }
    for( unsigned int imode = 1; imode < nmodes ; imode++ ) {
        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
        Jl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jl_[imode];
        Jr = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jr_[imode];
        Jt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Jt_[imode];
        Rho= ( static_cast<ElectroMagnAM *>( EMfields ) )->rho_AM_[imode];
        
        exp_mm_theta *= exp_m_theta ;
        //int zero_on_axis_l_rho  =  (imode != 0)*2;//=2 if El, Bl, Jl or rho is zero on axis and 0 if they are constant.
        int zero_on_axis_rt = (imode != 1)*2;     //=2 if Er, Et, Br, Bt, Jr, Jt is zero on axis and 0 if they are constant.
        
        *( ELoc+0*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[2], El, ip_, jp_ ) * exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Er, ip_, jp_ ) * exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Et, ip_, jp_ ) * exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[2], Bl, ip_, jp_ ) * exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Br, ip_, jp_ ) * exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Bt, ip_, jp_ ) * exp_mm_theta ) ;
        JLoc->x += std::real( compute( &coeffxp_[0], &coeffyp_[2              ], Jl, ip_, jp_ ) * exp_mm_theta ) ;
        JLoc->y += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Jr, ip_, jp_ ) * exp_mm_theta ) ;
        JLoc->z += std::real( compute( &coeffxp_[0], &coeffyp_[zero_on_axis_rt], Jt, ip_, jp_ ) * exp_mm_theta ) ;
        ( *RhoLoc ) += std::real( compute( &coeffxp_[0], &coeffyp_[2], Rho, ip_, jp_ )* exp_mm_theta ) ;
    }
    double delta2 = std::real( exp_m_theta ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta ) * *( ELoc+1*nparts ) + std::real( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta ) *  *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta ) * *( BLoc+1*nparts ) + std::real( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta ) * JLoc->y + std::imag( exp_m_theta ) * JLoc->z;
    JLoc->z = -std::imag( exp_m_theta ) * JLoc->y + std::real( exp_m_theta ) * JLoc->z;
    JLoc->y = delta2 ;
    
}

// Interpolator on another field than the basic ones
void InterpolatorAM1Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *Jxloc, double *Jyloc, double *Jzloc, double *Rholoc )
{
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*dl_inv_;
        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        complex<double> exp_m_theta, exp_mm_theta = 1. ;
        double rpn = r * dr_inv_ - 0.5;
        coeffs( xpn, rpn);
        if (r > 0){ 
            exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        } else {
            exp_m_theta = 1. ;
        }

        double Jx_ = 0., Jy_ = 0., Jz_ = 0., Rho_ = 0.;
        for( unsigned int imode = 0; imode < nmodes ; imode++ ) {
            int icoeff;
            if (imode == 0){ // if rho or longitudinal field
                icoeff = 0; // constant on axis
            } else {
                icoeff = 2; // zero on axis
            }
            cField2D *Jl  = static_cast<cField2D *>( *(field+4*imode+0) );
            cField2D *Jr  = static_cast<cField2D *>( *(field+4*imode+1) );
            cField2D *Jt  = static_cast<cField2D *>( *(field+4*imode+2) );
            cField2D *Rho = static_cast<cField2D *>( *(field+4*imode+3) );
            Jx_  += std::real( compute( &coeffxp_[1], &coeffyp_[icoeff], Jl , ip_, jp_ ) * exp_mm_theta );
            Jy_  += std::real( compute( &coeffxp_[1], &coeffyp_[(imode != 1)*2], Jr , ip_, jp_ ) * exp_mm_theta );
            Jz_  += std::real( compute( &coeffxp_[1], &coeffyp_[(imode != 1)*2], Jt , ip_, jp_ ) * exp_mm_theta );
            Rho_ += std::real( compute( &coeffxp_[1], &coeffyp_[icoeff], Rho, ip_, jp_ ) * exp_mm_theta );
            exp_mm_theta *= exp_m_theta;
        }
        Jxloc [ipart] = Jx_;
        Jyloc [ipart] = std::real( exp_m_theta ) * Jy_ + std::imag( exp_m_theta ) * Jz_;
        Jzloc [ipart] = -std::imag( exp_m_theta ) * Jy_ + std::real( exp_m_theta ) * Jz_;
        Rholoc[ipart] = Rho_;

    }

}

void InterpolatorAM1Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<double> *theta_old = &( smpi->dynamics_thetaold[ithread] );
    
    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );
        //Buffering of iol and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltar;
        ( *theta_old )[ipart] = atan2( particles.position( 2, ipart ), particles.position( 1, ipart ));
    }
    

}


// Interpolator specific to tracked particles. A selection of particles may be provided
void InterpolatorAM1Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
{
    ERROR( "To Do" );
}
