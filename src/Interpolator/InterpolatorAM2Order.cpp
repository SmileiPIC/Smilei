#include "InterpolatorAM2Order.h"

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
// Creator for InterpolatorAM2Order
// ---------------------------------------------------------------------------------------------------------------------
InterpolatorAM2Order::InterpolatorAM2Order( Params &params, Patch *patch ) : InterpolatorAM( patch )
{

    D_inv_[0] = 1.0/params.cell_length[0];
    D_inv_[1] = 1.0/params.cell_length[1];
    nmodes_ = params.nmodes;

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void InterpolatorAM2Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    // Calculate coeffs
    coeffs( xpn, rpn );

    //std::cout << "thetas = " << exp_m_theta_ << " " << 2.*std::real(exp_m_theta_) - exp_m_theta_ << std::endl;

    //Here we assume that mode 0 is real !!
    // Interpolation of El^(d,p)
    *( ELoc+0*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ ) );
    // Interpolation of Er^(p,d)
    *( ELoc+1*nparts ) = std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ ) );
    // Interpolation of Et^(p,p)
    *( ELoc+2*nparts ) = std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ ) );
    // Interpolation of Bl^(p,d)
    *( BLoc+0*nparts ) = std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ ) );
    // Interpolation of Br^(d,p)
    *( BLoc+1*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ ) );
    // Interpolation of Bt^(d,d)
    *( BLoc+2*nparts ) = std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ ) );
    //std::cout << "mode 0 " <<*( ELoc+0*nparts ) << " " <<*( ELoc+1*nparts )   << " " <<*( ELoc+2*nparts )<<" " <<*( BLoc+0*nparts ) << " " <<*( BLoc+1*nparts )   << " " <<*( BLoc+2*nparts ) << endl;

    for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];

        exp_mm_theta *= exp_m_theta_ ;

        *( ELoc+0*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ )* exp_mm_theta ) ;
        //std::cout << " compute Er mode 1 " << endl;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ )* exp_mm_theta ) ;
        //std::cout << " done compute Er mode 1 by adding " << std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ )* exp_mm_theta ) << " expmmtheta = " << exp_mm_theta << endl;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ )* exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ )* exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ )* exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ )* exp_mm_theta ) ;
    }
    //std::cout << "mode 1 " <<*( ELoc+0*nparts ) << " " <<*( ELoc+1*nparts )   << " " <<*( ELoc+2*nparts )<<" " <<*( BLoc+0*nparts ) << " " <<*( BLoc+1*nparts )   << " " <<*( BLoc+2*nparts ) << endl;

    //Translate field into the cartesian y,z coordinates
    double delta2 = std::real( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( ELoc+1*nparts ) + std::real( exp_m_theta_ ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta_ ) * *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta_ ) * *( BLoc+1*nparts ) + std::real( exp_m_theta_ ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;
    //std::cout << "standard final " <<  *( ELoc ) << " " <<  *( ELoc+nparts ) << " " <<  *( ELoc+2*nparts )<< " " <<  *( BLoc ) << " " <<  *( BLoc+nparts ) << " " <<  *( BLoc+2*nparts ) << std::endl;

} // END InterpolatorAM2Order

void InterpolatorAM2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    double xpn = particles.position( 0, ipart ) * D_inv_[0];
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * D_inv_[1];
    complex<double> exp_mm_theta = 1. ;

    // Calculate coeffs
    coeffs( xpn, rpn );

    int nparts( particles.numberOfParticles() );

    // Interpolation of El^(d,p)
    *( ELoc+0*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ ) );
    // Interpolation of Er^(p,d)
    *( ELoc+1*nparts ) = std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ ) );
    // Interpolation of Et^(p,p)
    *( ELoc+2*nparts ) = std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ ) );
    // Interpolation of Bl^(p,d)
    *( BLoc+0*nparts ) = std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ ) );
    // Interpolation of Br^(d,p)
    *( BLoc+1*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ ) );
    // Interpolation of Bt^(d,d)
    *( BLoc+2*nparts ) = std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ ) );
    // Interpolation of Jl^(d,p,p)
    JLoc->x = std::real( compute( &coeffxd_[1], &coeffyp_[1], Jl, id_, jp_ ) );
    // Interpolation of Jr^(p,d,p)
    JLoc->y = std::real( compute( &coeffxp_[1], &coeffyd_[1], Jr, ip_, jd_ ) );
    // Interpolation of Jt^(p,p,d)
    JLoc->z = std::real( compute( &coeffxp_[1], &coeffyp_[1], Jt, ip_, jp_ ) );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = std::real( compute( &coeffxp_[1], &coeffyp_[1], Rho, ip_, jp_ ) );

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

        exp_mm_theta *= exp_m_theta_ ;

        *( ELoc+0*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ ) * exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ ) * exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ ) * exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ ) * exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ ) * exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ ) * exp_mm_theta ) ;
        JLoc->x += std::real( compute( &coeffxd_[1], &coeffyp_[1], Jl, id_, jp_ ) * exp_mm_theta ) ;
        JLoc->y += std::real( compute( &coeffxp_[1], &coeffyd_[1], Jr, ip_, jd_ ) * exp_mm_theta ) ;
        JLoc->z += std::real( compute( &coeffxp_[1], &coeffyp_[1], Jt, ip_, jp_ ) * exp_mm_theta ) ;
        ( *RhoLoc ) += std::real( compute( &coeffxp_[1], &coeffyp_[1], Rho, ip_, jp_ )* exp_mm_theta ) ;
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

}

// Interpolator on another field than the basic ones
void InterpolatorAM2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *Jxloc, double *Jyloc, double *Jzloc, double *Rholoc )
{

    // **field points to the first field of the species of interest in EM->allFields
    // They are ordered as Jx0, Jy0, Jz0, Rho0, Jx1, Jy1, Jz1, Rho1, etc.


    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*D_inv_[0];
        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        double rpn = r * D_inv_[1];
        coeffs( xpn, rpn);
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
            Jx_  += std::real( compute( &coeffxd_[1], &coeffyp_[1], Jl , id_, jp_ ) * exp_mm_theta );
            Jy_  += std::real( compute( &coeffxp_[1], &coeffyd_[1], Jr , ip_, jd_ ) * exp_mm_theta );
            Jz_  += std::real( compute( &coeffxp_[1], &coeffyp_[1], Jt , ip_, jp_ ) * exp_mm_theta );
            Rho_ += std::real( compute( &coeffxp_[1], &coeffyp_[1], Rho, ip_, jp_ ) * exp_mm_theta );

            exp_mm_theta *= exp_m_theta_;
        }
        Jxloc [ipart] = Jx_;
        Jyloc [ipart] = std::real( exp_m_theta_ ) * Jy_ + std::imag( exp_m_theta_ ) * Jz_;
        Jzloc [ipart] = -std::imag( exp_m_theta_ ) * Jy_ + std::real( exp_m_theta_ ) * Jz_;
        Rholoc[ipart] = Rho_;
    }
}

void InterpolatorAM2Order::fieldsWrapper( ElectroMagn *EMfields,
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
        double coeffxp[3], coeffyp[3];
        double coeffxd[3], coeffyd[3];

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
        *( Epart+0*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyp[1], El, idx_d[0], idx_p[1] ) );
        // Interpolation of Er^(p,d)
        *( Epart+1*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyd[1], Er, idx_p[0], idx_d[1] ) );
        // Interpolation of Et^(p,p)
        *( Epart+2*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyp[1], Et, idx_p[0], idx_p[1] ) );
        // Interpolation of Bl^(p,d)
        *( Bpart+0*nparts+ipart ) = std::real( compute( &coeffxp[1], &coeffyd[1], Bl, idx_p[0], idx_d[1] ) );
        // Interpolation of Br^(d,p)
        *( Bpart+1*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyp[1], Br, idx_d[0], idx_p[1] ) );
        // Interpolation of Bt^(d,d)
        *( Bpart+2*nparts+ipart ) = std::real( compute( &coeffxd[1], &coeffyd[1], Bt, idx_d[0], idx_d[1] ) );

        for( unsigned int imode = 1; imode < nmodes_ ; imode++ ) {
            El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
            Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
            Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
            Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
            Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
            Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];

            exp_mm_theta_local *= exp_m_theta_local ;

            *( Epart+0*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyp[1], El, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
            *( Epart+1*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyd[1], Er, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
            *( Epart+2*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyp[1], Et, idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
            *( Bpart+0*nparts+ipart ) += std::real( compute( &coeffxp[1], &coeffyd[1], Bl, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
            *( Bpart+1*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyp[1], Br, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
            *( Bpart+2*nparts+ipart ) += std::real( compute( &coeffxd[1], &coeffyd[1], Bt, idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
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

    }

}

// Interpolator specific to tracked particles. A selection of particles may be provided
void InterpolatorAM2Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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


void InterpolatorAM2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
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

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        int idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        double coeffxd[3], coeffyd[3];
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
        ( *Epart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[1], El, idx_d[0], idx_p[1] ) );
        // Interpolation of Er^(p,d)
        ( *Epart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[1], Er, idx_p[0], idx_d[1] ) );
        // Interpolation of Et^(p,p)
        ( *Epart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyp[1], Et, idx_p[0], idx_p[1] ) );
        // Interpolation of Bl^(p,d)
        ( *Bpart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxp[1], &coeffyd[1], Bl, idx_p[0], idx_d[1] ) );
        // Interpolation of Br^(d,p)
        ( *Bpart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyp[1], Br, idx_d[0], idx_p[1] ) );
        // Interpolation of Bt^(d,d)
        ( *Bpart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxd[1], &coeffyd[1], Bt, idx_d[0], idx_d[1] ) );
        // Interpolation of Phi^(p,p)
        ( *PHIpart ) [ 0*nparts+ipart ]     = compute( &coeffxp[1], &coeffyp[1], Phi, idx_p[0], idx_p[1] ) ;
        // Interpolation of GradPhil^(p,p)
        ( *GradPHIpart ) [ 0*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[1], GradPhil, idx_p[0], idx_p[1] ) ;
        // Interpolation of GradPhir^(p,p)
        ( *GradPHIpart ) [ 1*nparts+ipart ] = compute( &coeffxp[1], &coeffyp[1], GradPhir, idx_p[0], idx_p[1] ) ;
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

            ( *Epart ) [ 0*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyp[1], El, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
            ( *Epart ) [ 1*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyd[1], Er, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
            ( *Epart ) [ 2*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyp[1], Et, idx_p[0], idx_p[1] )* exp_mm_theta_local ) ;
            ( *Bpart ) [ 0*nparts+ipart ] += std::real( compute( &coeffxp[1], &coeffyd[1], Bl, idx_p[0], idx_d[1] )* exp_mm_theta_local ) ;
            ( *Bpart ) [ 1*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyp[1], Br, idx_d[0], idx_p[1] )* exp_mm_theta_local ) ;
            ( *Bpart ) [ 2*nparts+ipart ] += std::real( compute( &coeffxd[1], &coeffyd[1], Bt, idx_d[0], idx_d[1] )* exp_mm_theta_local ) ;
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
        ( *eitheta_old)[ipart] =  2.*std::real(exp_m_theta_local) - exp_m_theta_local ;  //exp(i theta)

    }

} // END InterpolatorAM2Order

void InterpolatorAM2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
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

        int idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        double coeffxd[3], coeffyd[3];

        // Compute coefficients
        coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

        // only mode 0 is used

        // -------------------------
        // Interpolation of Phi_m^(p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], Phi_m2Dcyl, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhi_m^(p,p), l component
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhil_m2Dcyl, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhi_m^(p,p), r component
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhir_m2Dcyl, idx_p[0], idx_p[1] );

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
        ( *eitheta_old)[ipart] =  2.*std::real(exp_m_theta_local) - exp_m_theta_local ;  //exp(i theta)

    }

} // END InterpolatorAM2Order


void InterpolatorAM2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
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
    coeffs( xpn, rpn );

    // -------------------------
    // Interpolation of Env_A_abs_^(p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_A_abs_2Dcyl, ip_, jp_);

    // -------------------------
    // Interpolation of Env_Chi_^(p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_Chi_2Dcyl, ip_, jp_);

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_E_abs_2Dcyl, ip_, jp_);

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], Env_Ex_abs_2Dcyl, ip_, jp_);

} // END InterpolatorAM2Order

void InterpolatorAM2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
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

        int idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        double coeffxd[3], coeffyd[3];

        // Compute coefficients
        coeffs( xpn, rpn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

        // only mode 0 is used

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], EnvEabs, idx_p[0], idx_p[1] );

        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], EnvExabs, idx_p[0], idx_p[1] );
    }

} // END InterpolatorAM2Order
