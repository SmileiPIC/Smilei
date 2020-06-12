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
InterpolatorAM2Order::InterpolatorAM2Order( Params &params, Patch *patch ) : InterpolatorAM( params, patch )
{

    dl_inv_ = 1.0/params.cell_length[0];
    dr_inv_ = 1.0/params.cell_length[1];
    nmodes = params.nmodes;
    dr =  params.cell_length[1];
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
    double xpn = particles.position( 0, ipart ) * dl_inv_;
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * dr_inv_;
    exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ; //exp(-i theta)
    complex<double> exp_mm_theta = 1. ;                                                          //exp(-i m theta)
    // Calculate coeffs
    coeffs( xpn, rpn );
    
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
    
    for( unsigned int imode = 1; imode < nmodes ; imode++ ) {
        El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[imode];
        Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[imode];
        Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[imode];
        Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[imode];
        Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[imode];
        Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[imode];
        
        exp_mm_theta *= exp_m_theta ;
        
        *( ELoc+0*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ )* exp_mm_theta ) ;
        *( ELoc+1*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ )* exp_mm_theta ) ;
        *( ELoc+2*nparts ) += std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ )* exp_mm_theta ) ;
        *( BLoc+0*nparts ) += std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ )* exp_mm_theta ) ;
        *( BLoc+1*nparts ) += std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ )* exp_mm_theta ) ;
        *( BLoc+2*nparts ) += std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ )* exp_mm_theta ) ;
    }
    
    //Translate field into the cartesian y,z coordinates
    double delta2 = std::real( exp_m_theta ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta ) * *( ELoc+1*nparts ) + std::real( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta ) * *( BLoc+1*nparts ) + std::real( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;
    
} // END InterpolatorAM2Order

void InterpolatorAM2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    double rpn = r * dr_inv_;
    complex<double> exp_mm_theta = 1. ;
    
    // Calculate coeffs
    coeffs( xpn, rpn );
    
    int nparts( particles.size() );
    
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
void InterpolatorAM2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *Jxloc, double *Jyloc, double *Jzloc, double *Rholoc )
{

    // **field points to the first field of the species of interest in EM->allFields
    // They are ordered as Jx0, Jy0, Jz0, Rho0, Jx1, Jy1, Jz1, Rho1, etc.

    
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*dl_inv_;
        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        double rpn = r * dr_inv_;
        coeffs( xpn, rpn);
        complex<double> exp_m_theta = 1., exp_mm_theta = 1. ;
        if (r > 0) {
            exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        }



        double Jx_ = 0., Jy_ = 0., Jz_ = 0., Rho_ = 0.;
        for( unsigned int imode = 0; imode < nmodes ; imode++ ) {
            cField2D *Jl  = static_cast<cField2D *>( *(field+4*imode+0) );
            cField2D *Jr  = static_cast<cField2D *>( *(field+4*imode+1) );
            cField2D *Jt  = static_cast<cField2D *>( *(field+4*imode+2) );
            cField2D *Rho = static_cast<cField2D *>( *(field+4*imode+3) );
            Jx_  += std::real( compute( &coeffxd_[1], &coeffyp_[1], Jl , id_, jp_ ) * exp_mm_theta );
            Jy_  += std::real( compute( &coeffxp_[1], &coeffyd_[1], Jr , ip_, jd_ ) * exp_mm_theta );
            Jz_  += std::real( compute( &coeffxp_[1], &coeffyp_[1], Jt , ip_, jp_ ) * exp_mm_theta );
            Rho_ += std::real( compute( &coeffxp_[1], &coeffyp_[1], Rho, ip_, jp_ ) * exp_mm_theta );

            exp_mm_theta *= exp_m_theta;
        }
        Jxloc [ipart] = Jx_;
        Jyloc [ipart] = std::real( exp_m_theta ) * Jy_ + std::imag( exp_m_theta ) * Jz_;
        Jzloc [ipart] = -std::imag( exp_m_theta ) * Jy_ + std::real( exp_m_theta ) * Jz_;
        Rholoc[ipart] = Rho_;
    }
}

void InterpolatorAM2Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
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
void InterpolatorAM2Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
{
    if( selection ) {
    
        int nsel_tot = selection->size();
        for( int isel=0 ; isel<nsel_tot; isel++ ) {
            fields( EMfields, particles, ( *selection )[isel], offset, buffer+isel, buffer+isel+3*offset );
        }
        
    } else {
    
        int npart_tot = particles.size();
        for( int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields( EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset );
        }
    }
}


void InterpolatorAM2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    
    
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    
    // Interpolate E, B
    cField2D *El = ( static_cast<ElectroMagnAM *>( EMfields ) )->El_[0];
    cField2D *Er = ( static_cast<ElectroMagnAM *>( EMfields ) )->Er_[0];
    cField2D *Et = ( static_cast<ElectroMagnAM *>( EMfields ) )->Et_[0];
    cField2D *Bl = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bl_m[0];
    cField2D *Br = ( static_cast<ElectroMagnAM *>( EMfields ) )->Br_m[0];
    cField2D *Bt = ( static_cast<ElectroMagnAM *>( EMfields ) )->Bt_m[0];

    // Static cast of the envelope fields
    Field2D *Phi = static_cast<Field2D *>( EMfields->envelope->Phi_ );
    Field2D *GradPhil = static_cast<Field2D *>( EMfields->envelope->GradPhil_ );
    Field2D *GradPhir = static_cast<Field2D *>( EMfields->envelope->GradPhir_ );
    
    // auxiliary quantities    
    double delta2, xpn, r, rpn;
    int nparts = particles.size() ;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        xpn = particles.position( 0, ipart ) * dl_inv_;
        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        rpn = r * dr_inv_;
    
        // Calculate coeffs
        coeffs( xpn, rpn );
    
        
        // only mode 0 is used

        // Interpolation of El^(d,p)
        ( *Epart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxd_[1], &coeffyp_[1], El, id_, jp_ ) );
        // Interpolation of Er^(p,d)
        ( *Epart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxp_[1], &coeffyd_[1], Er, ip_, jd_ ) );
        // Interpolation of Et^(p,p)
        ( *Epart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxp_[1], &coeffyp_[1], Et, ip_, jp_ ) );
        // Interpolation of Bl^(p,d)
        ( *Bpart ) [ 0*nparts+ipart ]       = std::real( compute( &coeffxp_[1], &coeffyd_[1], Bl, ip_, jd_ ) );
        // Interpolation of Br^(d,p)
        ( *Bpart ) [ 1*nparts+ipart ]       = std::real( compute( &coeffxd_[1], &coeffyp_[1], Br, id_, jp_ ) );
        // Interpolation of Bt^(d,d)
        ( *Bpart ) [ 2*nparts+ipart ]       = std::real( compute( &coeffxd_[1], &coeffyd_[1], Bt, id_, jd_ ) );
        // Interpolation of Phi^(p,p)
        ( *PHIpart ) [ 0*nparts+ipart ]     = compute( &coeffxp_[1], &coeffyp_[1], Phi, ip_, jp_ ) ;
        // Interpolation of GradPhil^(p,p)
        ( *GradPHIpart ) [ 0*nparts+ipart ] = compute( &coeffxp_[1], &coeffyp_[1], GradPhil, ip_, jp_ ) ;
        // Interpolation of GradPhir^(p,p)
        ( *GradPHIpart ) [ 1*nparts+ipart ] = compute( &coeffxp_[1], &coeffyp_[1], GradPhir, ip_, jp_ ) ;
        // GradPhit = 0 in cylindrical symmetry
        ( *GradPHIpart ) [ 2*nparts+ipart ] = 0.;
   
        if (r > 0){ 
            exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        } else {
            exp_m_theta = 1. ;
        }

        // project on x,y,z, remember that GradPhit = 0 in cylindrical symmetry
        delta2 = std::real( exp_m_theta ) * ( *Epart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta ) * ( *Epart ) [ 2*nparts+ipart ];
        ( *Epart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta ) * ( *Epart ) [ 1*nparts+ipart ] + std::real( exp_m_theta ) * ( *Epart ) [ 2*nparts+ipart ];
        ( *Epart ) [ 1*nparts+ipart ] = delta2 ;
        delta2 = std::real( exp_m_theta ) * ( *Bpart ) [ 1*nparts+ipart ] + std::imag( exp_m_theta ) *  ( *Bpart ) [ 2*nparts+ipart ];
        ( *Bpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta ) * ( *Bpart ) [ 1*nparts+ipart ] + std::real( exp_m_theta ) * ( *Bpart ) [ 2*nparts+ipart ];
        ( *Bpart ) [ 1*nparts+ipart ] = delta2 ;

        delta2 = std::real( exp_m_theta ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ; 
        ( *GradPHIpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta ) * ( *GradPHIpart ) [ 1*nparts+ipart ] ;
        ( *GradPHIpart ) [ 1*nparts+ipart ] = delta2 ;

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltar;
        

    }
    
} // END InterpolatorAM2Order


void InterpolatorAM2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
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
    
    double r, delta2, xpn, rpn;
    //Loop on bin particles
    int nparts =  particles.size() ;
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
    
        // Normalized particle position
        xpn = particles.position( 0, ipart ) * dl_inv_;
        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        rpn = r * dr_inv_;
    
        // Compute coefficients
        coeffs( xpn, rpn );

        // only mode 0 is used
    
        // -------------------------
        // Interpolation of Phi_m^(p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], Phi_m2Dcyl, ip_, jp_ );
    
        // -------------------------
        // Interpolation of GradPhi_m^(p,p), l component
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhil_m2Dcyl, ip_, jp_ );
    
        // -------------------------
        // Interpolation of GradPhi_m^(p,p), r component
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], GradPhir_m2Dcyl, ip_, jp_ );
    
        // -------------------------
        // Interpolation of GradPhi_m^(p,p), theta component
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = 0.; // zero with cylindrical symmetry


        if (r > 0){ 
            exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        } else {
            exp_m_theta = 1. ;
        }


        // project on x,y,z, remember that GradPhit = 0 in cylindrical symmetry
        delta2 = std::real( exp_m_theta ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ; 
        ( *GradPHI_mpart ) [ 2*nparts+ipart ] = -std::imag( exp_m_theta ) * ( *GradPHI_mpart ) [ 1*nparts+ipart ] ;
        ( *GradPHI_mpart ) [ 1*nparts+ipart ] = delta2 ;

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
      
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltar;
      
    
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
    double xpn = particles.position( 0, ipart ) * dl_inv_;
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * dr_inv_;

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

void InterpolatorAM2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
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
        xpn = particles.position( 0, ipart ) * dl_inv_;
        r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        rpn = r * dr_inv_;
                                     
        // Compute coefficients
        coeffs( xpn, rpn );
 
        // only mode 0 is used
    
        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], EnvEabs, ip_, jp_ );
  
        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], EnvExabs, ip_, jp_ );
    
    }
    
    
} // END InterpolatorAM2Order




