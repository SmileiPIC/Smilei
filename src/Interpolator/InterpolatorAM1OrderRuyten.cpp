#include "InterpolatorAM1OrderRuyten.h"

#include <cmath>
#include <iostream>
#include <math.h>
#include "ElectroMagn.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include "Particles.h"
#include "LaserEnvelope.h"
#include <complex>
#include "dcomplex.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for InterpolatorAM1OrderRuyten
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorAM1OrderRuyten::InterpolatorAM1OrderRuyten( Params &params, Patch *patch ) : InterpolatorAM( patch )
{

    D_inv_[0] = 1.0/params.cell_length[0];
    D_inv_[1] = 1.0/params.cell_length[1];
    nmodes_ = params.nmodes;

}

// ---------------------------------------------------------------------------------------------------------------------
// 1st Order Interpolation of the fields with Ruyten correction at a the particle position (2 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void InterpolatorAM1OrderRuyten::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    double xpn = particles.position( 0, ipart ) * D_inv_[0];
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * D_inv_[1];
    exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ; //exp(-i theta)
    complex<double> exp_mm_theta = 1. ;                                                          //exp(-i m theta)
    // Compute coeffs
    int idx_p[2], idx_d[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[2];
    double coeffxd[3], coeffyd[2];
    coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

    //Here we assume that mode 0 is real !!
    // Interpolation of El^(d,p)
    *( ELoc+0*nparts )      = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
    // Interpolation of Er^(p,d)
    *( ELoc+1*nparts )      = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
    // Interpolation of Et^(p,p)
    *( ELoc+2*nparts )      = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
    // Interpolation of Bl^(p,d)
    *( BLoc+0*nparts )      = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
    // Interpolation of Br^(d,p)
    *( BLoc+1*nparts )      = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
    // Interpolation of Bt^(d,d)
    *( BLoc+2*nparts )      = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );

    for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];

        exp_mm_theta *= exp_m_theta_ ;

        *( ELoc+0*nparts ) += std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] )* exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] )* exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] )* exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] )* exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] )* exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] )* exp_mm_theta ) ;
    }

    //Translate field into the cartesian y,z coordinates
    double delta2 = std::real( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::real( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::real( exp_m_theta_ ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;

} // END InterpolatorAM1OrderRuyten

void InterpolatorAM1OrderRuyten::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );
    double *BLocyBTIS3;
    double *BLoczBTIS3;
    if (smpi->use_BTIS3){
        BLocyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][ipart] );
        BLoczBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][ipart] );
    }

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
    cField2D *Br_BTIS3;
    cField2D *Bt_BTIS3;
    if (smpi->use_BTIS3){
        Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[0];
        Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[0];
    }

    // Normalized particle position
    double xpn = particles.position( 0, ipart ) * D_inv_[0];
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * D_inv_[1];
    complex<double> exp_mm_theta = 1. ;

    // Calculate coeffs
    int idx_p[2], idx_d[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[2];
    double coeffxd[3], coeffyd[2];
    coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

    int nparts( particles.numberOfParticles() );

    *( ELoc+0*nparts )      = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
    // Interpolation of Er^(p,d)
    *( ELoc+1*nparts )      = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
    // Interpolation of Et^(p,p)
    *( ELoc+2*nparts )      = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
    // Interpolation of Bl^(p,d)
    *( BLoc+0*nparts )      = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
    // Interpolation of Br^(d,p)
    *( BLoc+1*nparts )      = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
    // Interpolation of Bt^(d,d)
    *( BLoc+2*nparts )      = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );

    // Interpolation of Jl^(d,p,p)
    JLoc->x = std::real( compute( &coeffxd[1], &coeffyp[0], Jl, idx_d[0], idx_p[1] ) );
    // Interpolation of Jr^(p,d,p)
    JLoc->y = std::real( compute( &coeffxp[1], &coeffyd[0], Jr, idx_p[0], idx_d[1]) );
    // Interpolation of Jt^(p,p,d)
    JLoc->z = std::real( compute( &coeffxp[1], &coeffyp[0], Jt, idx_p[0], idx_p[1]) );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = std::real( compute( &coeffxp[1], &coeffyp[0], Rho, idx_p[0], idx_p[1] ) );
    
    if (smpi->use_BTIS3){
        // BTIS fields, in the x direction they are centered as Ey and Ez
        // Interpolation of Br^(p,p) for BTIS
        *( BLocyBTIS3+0*nparts ) = std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] ) );
        // Interpolation of Bt^(p,d) for BTIS
        *( BLoczBTIS3+0*nparts ) = std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] ) );
    }

    if (r > 0){
        exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
    } else {
        exp_m_theta_ = 1. ;
    }
    for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
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
        
        if (smpi->use_BTIS3){
            Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[imode];
            Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[imode];
        }

        exp_mm_theta *= exp_m_theta_ ;

        *( ELoc+0*nparts ) += std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] )* exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] )* exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] )* exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] )* exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] )* exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] )* exp_mm_theta ) ;
        JLoc->x += std::real( compute( &coeffxd[1], &coeffyp[0], Jl, idx_d[0], idx_p[1] ) * exp_mm_theta ) ;
        JLoc->y += std::real( compute( &coeffxp[1], &coeffyd[0], Jr, idx_p[0], idx_d[1] ) * exp_mm_theta ) ;
        JLoc->z += std::real( compute( &coeffxp[1], &coeffyp[0], Jt, idx_p[0], idx_p[1] ) * exp_mm_theta ) ;
        ( *RhoLoc ) += std::real( compute( &coeffxp[1], &coeffyp[0], Rho, idx_p[0], idx_p[1] )* exp_mm_theta ) ;
        
        if (smpi->use_BTIS3){
            // BTIS3 fields, in the x direction they are centered as Ey and Ez
            *( BLocyBTIS3+0*nparts ) += std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] )* exp_mm_theta );
            *( BLoczBTIS3+0*nparts ) += std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] )* exp_mm_theta );
        }
    }
    double delta2 = std::real( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::real( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta_ ) *  *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::real( exp_m_theta_ ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta_ ) * JLoc->y + std::imag( exp_m_theta_ ) * JLoc->z;
    JLoc->z = -std::imag( exp_m_theta_ ) * JLoc->y + std::real( exp_m_theta_ ) * JLoc->z;
    JLoc->y = delta2 ;
    if (smpi->use_BTIS3){
    delta2 = std::real( exp_m_theta_ ) * *( BLocyBTIS3+0*nparts ) + std::imag( exp_m_theta_ ) * *( BLoczBTIS3+0*nparts );
    *( BLoczBTIS3+0*nparts ) = -std::imag( exp_m_theta_ ) * *( BLocyBTIS3+0*nparts ) + std::real( exp_m_theta_ ) * *( BLoczBTIS3+0*nparts );
    *( BLocyBTIS3+0*nparts ) = delta2 ;
    }

}

// Interpolator on another field than the basic ones
void InterpolatorAM1OrderRuyten::oneField( Field **field, Particles &particles, int *istart, int *iend, double *Jxloc, double *Jyloc, double *Jzloc, double *Rholoc )
{

    // **field points to the first field of the species of interest in EM->allFields
    // They are ordered as Jx0, Jy0, Jz0, Rho0, Jx1, Jy1, Jz1, Rho1, etc.

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*D_inv_[0];
        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        double rpn = r * D_inv_[1];

        int idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[2];
        double coeffxd[3], coeffyd[2];
        coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

        complex<double> exp_m_theta_ = 1., exp_mm_theta = 1. ;
        if (r > 0) {
            exp_m_theta_ = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        }

        double Jx_ = 0., Jy_ = 0., Jz_ = 0., Rho_ = 0.;
        for( unsigned int imode = 0; imode < nmodes_ ; imode++ ) {
            cField2D *Jl  = static_cast<cField2D *>( *(field+4*imode+0) );
            cField2D *Jr  = static_cast<cField2D *>( *(field+4*imode+1) );
            cField2D *Jt  = static_cast<cField2D *>( *(field+4*imode+2) );
            cField2D *Rho = static_cast<cField2D *>( *(field+4*imode+3) );
            Jx_  += std::real( compute( &coeffxd[1], &coeffyp[0], Jl , idx_d[0], idx_p[1] ) * exp_mm_theta );
            Jy_  += std::real( compute( &coeffxp[1], &coeffyd[0], Jr , idx_p[0], idx_d[1]  ) * exp_mm_theta );
            Jz_  += std::real( compute( &coeffxp[1], &coeffyp[0], Jt , idx_p[0], idx_p[1] ) * exp_mm_theta );
            Rho_ += std::real( compute( &coeffxp[1], &coeffyp[0], Rho, idx_p[0], idx_p[1] ) * exp_mm_theta );

            exp_mm_theta *= exp_m_theta_;
        }
        Jxloc [ipart] = Jx_;
        Jyloc [ipart] = std::real( exp_m_theta_ ) * Jy_ + std::imag( exp_m_theta_ ) * Jz_;
        Jzloc [ipart] = -std::imag( exp_m_theta_ ) * Jy_ + std::real( exp_m_theta_ ) * Jz_;
        Rholoc[ipart] = Rho_;
    }
}

void InterpolatorAM1OrderRuyten::fieldsWrapper( ElectroMagn *EMfields,
                                          Particles &particles,
                                          SmileiMPI *smpi,
                                          int *istart,
                                          int *iend,
                                          int ithread,
                                          unsigned int,
                                          int )
{

    double *Epart = &( smpi->dynamics_Epart[ithread][0] );
    double *Bpart = &( smpi->dynamics_Bpart[ithread][0] );
    int *iold = &( smpi->dynamics_iold[ithread][0] );
    double *delta = &( smpi->dynamics_deltaold[ithread][0] );
    std::complex<double> *eitheta_old = &( smpi->dynamics_eithetaold[ithread][0] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    
    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            complex<double> exp_m_theta_local ;
            complex<double> exp_mm_theta_local = 1. ;

            // Normalized particle position
            double xpn = particles.position( 0, ipart ) * D_inv_[0];
            double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
            double rpn = r * D_inv_[1];
            exp_m_theta_local = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ; //exp(-i theta)
                                                                   //exp(-i m theta)

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[2];
            double coeffxd[3], coeffyd[2];
            coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // Static cast of the electromagnetic fields, mode 0
            cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
            cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
            cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
            cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
            cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
            cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];

            //Here we assume that mode 0 is real !!
            // Interpolation of El^(d,p)
            *( Epart+0*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
            // Interpolation of Er^(p,d)
            *( Epart+1*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
            // Interpolation of Et^(p,p)
            *( Epart+2*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
            // Interpolation of Bl^(p,d)
            *( Bpart+0*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
            // Interpolation of Br^(d,p)
            *( Bpart+1*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
            // Interpolation of Bt^(d,d)
            *( Bpart+2*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );

            for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
                El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
                Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
                Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
                Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
                Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
                Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];

                exp_mm_theta_local *= exp_m_theta_local ;

                *( Epart+0*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Epart+1*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                *( Epart+2*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Bpart+0*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                *( Bpart+1*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Bpart+2*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
            }

            //Translate field into the cartesian y,z coordinates
            double delta2 = std::real( exp_m_theta_local ) * *( Epart+1*nparts+ipart ) + std::imag( exp_m_theta_local ) * *( Epart+2*nparts+ipart );
            *( Epart+2*nparts+ipart ) = -std::imag( exp_m_theta_local ) * *( Epart+1*nparts+ipart ) + std::real( exp_m_theta_local ) * *( Epart+2*nparts+ipart );
            *( Epart+1*nparts+ipart ) = delta2 ;
            delta2 = std::real( exp_m_theta_local ) * *( Bpart+1*nparts+ipart ) + std::imag( exp_m_theta_local ) * *( Bpart+2*nparts+ipart );
            *( Bpart+2*nparts+ipart ) = -std::imag( exp_m_theta_local ) * *( Bpart+1*nparts+ipart ) + std::real( exp_m_theta_local ) * *( Bpart+2*nparts+ipart );
            *( Bpart+1*nparts+ipart ) = delta2 ;

            // store indices and delta
            *( iold+0*nparts+ipart)  = idx_p[0];
            *( iold+1*nparts+ipart)  = idx_p[1];
            *( delta+0*nparts+ipart) = delta_p[0];
            *( delta+1*nparts+ipart) = delta_p[1];
            *( eitheta_old+ipart )     = 2.*std::real(exp_m_theta_local) - exp_m_theta_local ;  //exp(i theta)

        } // end ipart loop
      
    } else { // with B-TIS3 interpolation
      
        double *BLocyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][0] );
        double *BLoczBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][0] );
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            complex<double> exp_m_theta_local ;
            complex<double> exp_mm_theta_local = 1. ;

            // Normalized particle position
            double xpn = particles.position( 0, ipart ) * D_inv_[0];
            double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
            double rpn = r * D_inv_[1];
            exp_m_theta_local = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ; //exp(-i theta)
                                                                   //exp(-i m theta)

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[2];
            double coeffxd[3], coeffyd[2];

            coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // Static cast of the electromagnetic fields, mode 0
            cField2D *El       = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
            cField2D *Er       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
            cField2D *Et       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
            cField2D *Bl       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
            cField2D *Br       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
            cField2D *Bt       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
            cField2D *Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[0];
            cField2D *Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[0];

            //Here we assume that mode 0 is real !!
            // Interpolation of El^(d,p)
            *( Epart+0*nparts+ipart )      = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
            // Interpolation of Er^(p,d)
            *( Epart+1*nparts+ipart )      = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
            // Interpolation of Et^(p,p)
            *( Epart+2*nparts+ipart )      = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
            // Interpolation of Bl^(p,d)
            *( Bpart+0*nparts+ipart )      = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
            // Interpolation of Br^(d,p)
            *( Bpart+1*nparts+ipart )      = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
            // Interpolation of Bt^(d,d)
            *( Bpart+2*nparts+ipart )      = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );
            
            // BTIS3 fields, in the x direction they are centered as Ey and Ez
            // Interpolation of Br^(p,p) for BTIS
            *( BLocyBTIS3+0*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] ) );
            // Interpolation of Bt^(p,d) for BTIS
            *( BLoczBTIS3+0*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] ) );

            for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
                El       = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
                Er       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
                Et       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
                Bl       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
                Br       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
                Bt       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
                Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[imode];
                Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[imode];

                exp_mm_theta_local *= exp_m_theta_local ;

                *( Epart+0*nparts+ipart )      += std::real( compute( &coeffxd[1], &coeffyp[0], El      , idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Epart+1*nparts+ipart )      += std::real( compute( &coeffxp[1], &coeffyd[0], Er      , idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                *( Epart+2*nparts+ipart )      += std::real( compute( &coeffxp[1], &coeffyp[0], Et      , idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Bpart+0*nparts+ipart )      += std::real( compute( &coeffxp[1], &coeffyd[0], Bl      , idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                *( Bpart+1*nparts+ipart )      += std::real( compute( &coeffxd[1], &coeffyp[0], Br      , idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                *( Bpart+2*nparts+ipart )      += std::real( compute( &coeffxd[1], &coeffyd[0], Bt      , idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
                // BTIS3 fields, in the x direction they are centered as Ey and Ez
                *( BLocyBTIS3+0*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] )* exp_mm_theta_local );
                *( BLoczBTIS3+0*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] )* exp_mm_theta_local );
            }

            //Translate field into the cartesian y,z coordinates
            double delta2 = std::real( exp_m_theta_local ) * *( Epart+1*nparts+ipart ) + std::imag( exp_m_theta_local ) * *( Epart+2*nparts+ipart );
            *( Epart+2*nparts+ipart ) = -std::imag( exp_m_theta_local ) * *( Epart+1*nparts+ipart ) + std::real( exp_m_theta_local ) * *( Epart+2*nparts+ipart );
            *( Epart+1*nparts+ipart ) = delta2 ;
            delta2 = std::real( exp_m_theta_local ) * *( Bpart+1*nparts+ipart ) + std::imag( exp_m_theta_local ) * *( Bpart+2*nparts+ipart );
            *( Bpart+2*nparts+ipart ) = -std::imag( exp_m_theta_local ) * *( Bpart+1*nparts+ipart ) + std::real( exp_m_theta_local ) * *( Bpart+2*nparts+ipart );
            *( Bpart+1*nparts+ipart ) = delta2 ;
            delta2 = std::real( exp_m_theta_local ) * *( BLocyBTIS3+0*nparts+ipart ) + std::imag( exp_m_theta_local ) * *( BLoczBTIS3+0*nparts+ipart );
            *( BLoczBTIS3+0*nparts+ipart ) = -std::imag( exp_m_theta_local ) * *( BLocyBTIS3+0*nparts+ipart ) + std::real( exp_m_theta_local ) * *( BLoczBTIS3+0*nparts+ipart );
            *( BLocyBTIS3+0*nparts+ipart ) = delta2 ;

            // store indices and delta
            *( iold+0*nparts+ipart)  = idx_p[0];
            *( iold+1*nparts+ipart)  = idx_p[1];
            *( delta+0*nparts+ipart) = delta_p[0];
            *( delta+1*nparts+ipart) = delta_p[1];
            *( eitheta_old+ipart )     = 2.*std::real(exp_m_theta_local) - exp_m_theta_local ;  //exp(i theta)

        } // end ipart loop
      
    } // end with B-TIS3 interpolation

}

// Interpolator specific to tracked particles. A selection of particles may be provided
void InterpolatorAM1OrderRuyten::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
{
    if( selection ) {

        int nsel_tot = selection->size();
        for( int isel=0 ; isel<nsel_tot; isel++ ) {
            fields( EMfields, particles, ( *selection )[isel], offset, buffer+isel, buffer+isel+3*offset );
        }

    } else {

        int npart_tot = particles.numberOfParticles();
        for( int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields( EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset );
        }
    }
}


void InterpolatorAM1OrderRuyten::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<std::complex<double>> *eitheta_old = &( smpi->dynamics_eithetaold[ithread] );

    // Static cast of the envelope fields
    Field2D *Phi = static_cast<Field2D *>( EMfields->envelope->Phi_ );
    Field2D *GradPhil = static_cast<Field2D *>( EMfields->envelope->GradPhil_ );
    Field2D *GradPhir = static_cast<Field2D *>( EMfields->envelope->GradPhir_ );

    // auxiliary quantities
    double delta2, xpn, r, rpn;
    int nparts = particles.numberOfParticles() ;

    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[2];
            double coeffxd[3], coeffyd[2];
            complex<double> exp_m_theta_local ;
            complex<double> exp_mm_theta_local = 1. ;

            // Normalized particle position
            xpn = particles.position( 0, ipart ) * D_inv_[0];
            r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
            rpn = r * D_inv_[1];

            // Compute coefficients
            coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // mode 0 is treated first

            // Interpolate E, B
            cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
            cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
            cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
            cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
            cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
            cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];

            // Interpolation of El^(d,p)
            ( *Epart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
            // Interpolation of Er^(p,d)
            ( *Epart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
            // Interpolation of Et^(p,p)
            ( *Epart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
            // Interpolation of Bl^(p,d)
            ( *Bpart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
            // Interpolation of Br^(d,p)
            ( *Bpart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
            // Interpolation of Bt^(d,d)
            ( *Bpart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );
            // Interpolation of Phi^(p,p)
            ( *PHIpart ) [ 0*nparts+ipart ]     = compute( &coeffxp[1], &coeffyp[0], Phi, idx_p[0], idx_p[1] ) ;
            // Interpolation of GradPhil^(p,p)
            ( *GradPHIpart ) [ 0*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[0], GradPhil, idx_p[0], idx_p[1] ) ;
            // Interpolation of GradPhir^(p,p)
            ( *GradPHIpart ) [ 1*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[0], GradPhir, idx_p[0], idx_p[1] ) ;
            // GradPhit = 0 in cylindrical symmetry
            ( *GradPHIpart ) [ 2*nparts+ipart ] = 0.;

            if (r > 0){
                exp_m_theta_local = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
            } else {
                exp_m_theta_local = 1. ;
            }
            for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
                El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
                Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
                Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
                Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
                Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
                Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];

                exp_mm_theta_local *= exp_m_theta_local ;
                
                ( *Epart ) [ 0*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Epart ) [ 1*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                ( *Epart ) [ 2*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 0*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 1*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 2*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
            }

            // project on x,y,z, remember that GradPhit = 0 in cylindrical symmetry
            delta2 = std::real( exp_m_theta_local ) * ( *Epart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_local ) * ( *Epart ) [ 2*nparts+ipart ];
            ( *Epart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *Epart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_local ) * ( *Epart ) [ 2*nparts+ipart ];
            ( *Epart ) [ 1*nparts+ipart ] = delta2 ;
            delta2 = std::real( exp_m_theta_local ) * ( *Bpart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_local ) *  ( *Bpart ) [ 2*nparts+ipart ];
            ( *Bpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *Bpart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_local ) * ( *Bpart ) [ 2*nparts+ipart ];
            ( *Bpart ) [ 1*nparts+ipart ] = delta2 ;

            delta2 = std::real( exp_m_theta_local ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
            ( *GradPHIpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
            ( *GradPHIpart ) [ 1*nparts+ipart ] = delta2 ;

            //Buffering of iold and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];
            ( *eitheta_old)[ipart] =  std::conj(exp_m_theta_local);  //exp(i theta)

        } // end ipart loop
      
    } else { // with B-TIS3 interpolation

        std::vector<double> *BLocyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread] );
        std::vector<double> *BLoczBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread] );
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[2];
            double coeffxd[3], coeffyd[2];
            complex<double> exp_m_theta_local ;
            complex<double> exp_mm_theta_local = 1. ;

            // Normalized particle position
            xpn = particles.position( 0, ipart ) * D_inv_[0];
            r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
            rpn = r * D_inv_[1];

            // Compute coefficients
            coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // mode 0 is treated first

            // Interpolate E, B
            cField2D *El       = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
            cField2D *Er       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
            cField2D *Et       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
            cField2D *Bl       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
            cField2D *Br       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
            cField2D *Bt       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];
            cField2D *Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[0];
            cField2D *Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[0];

            // Interpolation of El^(d,p)
            ( *Epart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[0], El, idx_d[0], idx_p[1] ) );
            // Interpolation of Er^(p,d)
            ( *Epart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[0], Er, idx_p[0], idx_d[1] ) );
            // Interpolation of Et^(p,p)
            ( *Epart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyp[0], Et, idx_p[0], idx_p[1] ) );
            // Interpolation of Bl^(p,d)
            ( *Bpart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[0], Bl, idx_p[0], idx_d[1] ) );
            // Interpolation of Br^(d,p)
            ( *Bpart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[0], Br, idx_d[0], idx_p[1] ) );
            // Interpolation of Bt^(d,d)
            ( *Bpart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyd[0], Bt, idx_d[0], idx_d[1] ) );
            // Interpolation of Phi^(p,p)
            ( *PHIpart ) [ 0*nparts+ipart ]     = compute( &coeffxp[1], &coeffyp[0], Phi, idx_p[0], idx_p[1] ) ;
            // Interpolation of GradPhil^(p,p)
            ( *GradPHIpart ) [ 0*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[0], GradPhil, idx_p[0], idx_p[1] ) ;
            // Interpolation of GradPhir^(p,p)
            ( *GradPHIpart ) [ 1*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[0], GradPhir, idx_p[0], idx_p[1] ) ;
            // GradPhit = 0 in cylindrical symmetry
            ( *GradPHIpart ) [ 2*nparts+ipart ] = 0.;
            // BTIS3 fields, in the x direction they are centered as Ey and Ez
            // Interpolation of Br^(p,p) for BTIS
            ( *BLocyBTIS3)[ 0*nparts+ipart ]    = std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] ) );
            // Interpolation of Bt^(p,d) for BTIS
            ( *BLoczBTIS3)[ 0*nparts+ipart ]    = std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] ) );

            if (r > 0){
                exp_m_theta_local = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
            } else {
                exp_m_theta_local = 1. ;
            }
            for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
                El       = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
                Er       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
                Et       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
                Bl       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
                Br       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
                Bt       = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
                Br_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_mBTIS3[imode];
                Bt_BTIS3 = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_mBTIS3[imode];

                exp_mm_theta_local *= exp_m_theta_local ;
              
                ( *Epart ) [ 0*nparts+ipart ]    += std::real( compute( &coeffxd[1], &coeffyp[0], El      , idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Epart ) [ 1*nparts+ipart ]    += std::real( compute( &coeffxp[1], &coeffyd[0], Er      , idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                ( *Epart ) [ 2*nparts+ipart ]    += std::real( compute( &coeffxp[1], &coeffyp[0], Et      , idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 0*nparts+ipart ]    += std::real( compute( &coeffxp[1], &coeffyd[0], Bl      , idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 1*nparts+ipart ]    += std::real( compute( &coeffxd[1], &coeffyp[0], Br      , idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
                ( *Bpart ) [ 2*nparts+ipart ]    += std::real( compute( &coeffxd[1], &coeffyd[0], Bt      , idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
                // BTIS3 fields, in the x direction they are centered as Ey and Ez
                ( *BLocyBTIS3)[ 0*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyp[0], Br_BTIS3, idx_p[0], idx_p[1] )* exp_mm_theta_local );
                ( *BLoczBTIS3)[ 0*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyd[0], Bt_BTIS3, idx_p[0], idx_d[1] )* exp_mm_theta_local );
            }

            // project on x,y,z, remember that GradPhit = 0 in cylindrical symmetry
            delta2 = std::real( exp_m_theta_local ) * ( *Epart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_local ) * ( *Epart ) [ 2*nparts+ipart ];
            ( *Epart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *Epart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_local ) * ( *Epart ) [ 2*nparts+ipart ];
            ( *Epart ) [ 1*nparts+ipart ] = delta2 ;
            delta2 = std::real( exp_m_theta_local ) * ( *Bpart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta_local ) *  ( *Bpart ) [ 2*nparts+ipart ];
            ( *Bpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *Bpart ) [ 1*nparts+ipart ] + std::real( exp_m_theta_local ) * ( *Bpart ) [ 2*nparts+ipart ];
            ( *Bpart ) [ 1*nparts+ipart ] = delta2 ;

            delta2 = std::real( exp_m_theta_local ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
            ( *GradPHIpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
            ( *GradPHIpart ) [ 1*nparts+ipart ] = delta2 ;
            
            delta2 = std::real( exp_m_theta_local ) * ( *BLocyBTIS3) [ 0*nparts+ipart ] + std::imag( exp_m_theta_local ) * ( *BLoczBTIS3) [ 0*nparts+ipart ];
            ( *BLoczBTIS3) [0*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *BLocyBTIS3)[ 0*nparts+ipart ] + std::real( exp_m_theta_local ) * ( *BLoczBTIS3) [ 0*nparts+ipart ];
            ( *BLocyBTIS3) [0*nparts+ipart ] = delta2 ;

            //Buffering of iold and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];
            ( *eitheta_old)[ipart] =  std::conj(exp_m_theta_local);  //exp(i theta)

        } // end ipart loop
      
    } // end with B-TIS3 interpolation
    

} // END InterpolatorAM1OrderRuyten

void InterpolatorAM1OrderRuyten::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // // Static cast of the envelope fields
    // Static cast of the envelope fields
    Field2D *Phi_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->Phi_m );
    Field2D *GradPhil_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->GradPhil_m );
    Field2D *GradPhir_m2Dcyl = static_cast<Field2D *>( EMfields->envelope->GradPhir_m );

    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    std::vector<std::complex<double>> *eitheta_old = &( smpi->dynamics_eithetaold[ithread] );

    double r, delta2, xpn, rpn;

    complex<double> exp_m_theta_local ;
    //Loop on bin particles
    int nparts =  particles.numberOfParticles() ;
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        xpn = particles.position( 0, ipart ) * D_inv_[0];
        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        rpn = r * D_inv_[1];

        int idx_p[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[2];

        // Compute coefficients
        coeffs( xpn, rpn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

        // only mode 0 is used

        // -------------------------
        // Interpolation of Phi_m^(p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp[1], &coeffyp[0], Phi_m2Dcyl, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhi_m^(p,p), l component
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[0], GradPhil_m2Dcyl, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhi_m^(p,p), r component
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[0], GradPhir_m2Dcyl, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhi_m^(p,p), theta component
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = 0.; // zero with cylindrical symmetry

        if (r > 0){
            exp_m_theta_local = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        } else {
            exp_m_theta_local = 1. ;
        }

        // project on x,y,z, remember that GradPhit = 0 in cylindrical symmetry
        delta2 = std::real( exp_m_theta_local ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ;
        ( *GradPHI_mpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta_local ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ;
        ( *GradPHI_mpart ) [ 1*nparts+ipart ] = delta2 ;

        //Buffering of iold and delta
        // store indices and delta
        ( *iold )[0*nparts+ipart]  = idx_p[0];
        ( *iold )[1*nparts+ipart]  = idx_p[1];
        ( *delta )[0*nparts+ipart] = delta_p[0];
        ( *delta )[1*nparts+ipart] = delta_p[1];
        ( *eitheta_old)[ipart] =  std::conj(exp_m_theta_local);  //exp(i theta)

    }

} // END InterpolatorAM1OrderRuyten


void InterpolatorAM1OrderRuyten::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    // Static cast of the electromagnetic fields
    Field2D *Env_A_abs_2Dcyl  = static_cast<Field2D *>( EMfields->Env_A_abs_ );
    Field2D *Env_Chi_2Dcyl    = static_cast<Field2D *>( EMfields->Env_Chi_ );
    Field2D *Env_E_abs_2Dcyl  = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    Field2D *Env_Ex_abs_2Dcyl = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart ) * D_inv_[0];
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * D_inv_[1];

    // Compute coefficients
    int idx_p[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[2];
    coeffs( xpn, rpn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

    // -------------------------
    // Interpolation of Env_A_abs_^(p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp[1], &coeffyp[0], Env_A_abs_2Dcyl, idx_p[0], idx_p[1]);

    // -------------------------
    // Interpolation of Env_Chi_^(p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp[1], &coeffyp[0], Env_Chi_2Dcyl, idx_p[0], idx_p[1]);

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp[1], &coeffyp[0], Env_E_abs_2Dcyl, idx_p[0], idx_p[1]);

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp[1], &coeffyp[0], Env_Ex_abs_2Dcyl, idx_p[0], idx_p[1]);

} // END InterpolatorAM1OrderRuyten

void InterpolatorAM1OrderRuyten::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field2D *EnvEabs  = static_cast<Field2D*>( EMfields->Env_E_abs_ );
    Field2D *EnvExabs = static_cast<Field2D*>( EMfields->Env_Ex_abs_ );

    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );

    double xpn,rpn,r;

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        xpn = particles.position( 0, ipart ) * D_inv_[0];
        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        rpn = r * D_inv_[1];

        // Compute coefficients
        int idx_p[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[2];
        coeffs( xpn, rpn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

        // only mode 0 is used

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[0], EnvEabs, idx_p[0], idx_p[1] );

        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[0], EnvExabs, idx_p[0], idx_p[1] );
    }

} // END InterpolatorAM1OrderRuyten
