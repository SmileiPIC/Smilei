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
void InterpolatorAM2Order::oneField( Field *field, Particles &particles, int *istart, int *iend, double *FieldLoc )
{
    cField2D *F = static_cast<cField2D *>( field );
    nmodes = std::stoi(&(F->name.back())); //Works only if number_of_AM < 10

    double *coeffx = field->isDual( 0 ) ? &coeffxd_[1] : &coeffxp_[1];
    double *coeffy = field->isDual( 1 ) ? &coeffyd_[1] : &coeffyp_[1];
    int *i = field->isDual( 0 ) ? &id_ : &ip_;
    int *j = field->isDual( 1 ) ? &jd_ : &jp_;
    
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        complex<double> exp_m_theta, exp_mm_theta = 1. ;
        double xpn = particles.position( 0, ipart )*dl_inv_;
        double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
        double rpn = r * dr_inv_;
        if (r > 0){ 
            exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
        } else {
            exp_m_theta = 1. ;
        }
        for (unsigned int imode=0; imode < nmodes; imode++){
            exp_mm_theta *= exp_m_theta;
        }
        coeffs( xpn, rpn);
        FieldLoc[ipart] = std::real( compute( coeffx, coeffy, F, *i, *j) * exp_mm_theta);
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
    ERROR( "To Do" );
}


void InterpolatorAM2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    int ipart = *istart;
    
    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    double *PHILoc        = &( smpi->dynamics_PHIpart[ithread][ipart] );
    double *GradPHILoc    = &( smpi->dynamics_GradPHIpart[ithread][ipart] );
    
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
    
    
    // Normalized particle position
    double xpn = particles.position( 0, ipart ) * dl_inv_;
    double r = sqrt( particles.position( 1, ipart )*particles.position( 1, ipart )+particles.position( 2, ipart )*particles.position( 2, ipart ) ) ;
    double rpn = r * dr_inv_;
    
    
    // Calculate coeffs
    coeffs( xpn, rpn );
    
    int nparts( particles.size() );
    

    // only mode 0 is used

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
    // Interpolation of Phi^(p,p)
    *( PHILoc+0*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], Phi, id_, jp_ ) );
    // Interpolation of GradPhil^(p,p)
    *( GradPHILoc+0*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], GradPhil, ip_, jp_ ) );
    // Interpolation of GradPhir^(p,p)
    *( GradPHILoc+1*nparts ) = std::real( compute( &coeffxd_[1], &coeffyp_[1], GradPhir, ip_, jp_ ) );
    // GradPhit = 0 in cylindrical symmetry
   
    if (r > 0){ 
        exp_m_theta = ( particles.position( 1, ipart ) - Icpx * particles.position( 2, ipart ) ) / r ;
    } else {
        exp_m_theta = 1. ;
    }

    
    // project on x,y,z
    double delta2 = std::real( exp_m_theta ) * *( ELoc+1*nparts ) + std::imag( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+2*nparts ) = -std::imag( exp_m_theta ) * *( ELoc+1*nparts ) + std::real( exp_m_theta ) * *( ELoc+2*nparts );
    *( ELoc+1*nparts ) = delta2 ;
    delta2 = std::real( exp_m_theta ) * *( BLoc+1*nparts ) + std::imag( exp_m_theta ) *  *( BLoc+2*nparts );
    *( BLoc+2*nparts ) = -std::imag( exp_m_theta ) * *( BLoc+1*nparts ) + std::real( exp_m_theta ) * *( BLoc+2*nparts );
    *( BLoc+1*nparts ) = delta2 ;

    delta2 = std::real( exp_m_theta ) * *( GradPHILoc+1*nparts ) ; 
    *( GradPHILoc+2*nparts ) = -std::imag( exp_m_theta ) * *( GradPHILoc+1*nparts ) ;
    *( GradPHILoc+1*nparts ) = delta2 ;
    
} // END Interpolator3D2Order


void InterpolatorAM2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // // Static cast of the envelope fields
    // Field3D *Phi_m3D = static_cast<Field3D *>( EMfields->envelope->Phi_m );
    // Field3D *GradPhix_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhix_m );
    // Field3D *GradPhiy_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhiy_m );
    // Field3D *GradPhiz_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhiz_m );
    // 
    // std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    // std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );
    // 
    // std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    // std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    // 
    // //Loop on bin particles
    // int nparts( particles.size() );
    // for( int ipart=*istart ; ipart<*iend; ipart++ ) {
    // 
    //     // Normalized particle position
    //     double xpn = particles.position( 0, ipart )*dx_inv_;
    //     double ypn = particles.position( 1, ipart )*dy_inv_;
    //     double zpn = particles.position( 2, ipart )*dz_inv_;
    // 
    //     coeffs( xpn, ypn, zpn );
    // 
    //     // -------------------------
    //     // Interpolation of Phi_m^(p,p,p)
    //     // -------------------------
    //     ( *PHI_mpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Phi_m3D, ip_, jp_, kp_ );
    // 
    //     // -------------------------
    //     // Interpolation of GradPhix_m^(p,p,p)
    //     // -------------------------
    //     ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhix_m3D, ip_, jp_, kp_ );
    // 
    //     // -------------------------
    //     // Interpolation of GradPhiy_m^(p,p,p)
    //     // -------------------------
    //     ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiy_m3D, ip_, jp_, kp_ );
    // 
    //     // -------------------------
    //     // Interpolation of GradPhiz_m^(p,p,p)
    //     // -------------------------
    //     ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiz_m3D, ip_, jp_, kp_ );
    // 
    //     //Buffering of iold and delta
    //     ( *iold )[ipart+0*nparts]  = ip_;
    //     ( *iold )[ipart+1*nparts]  = jp_;
    //     ( *iold )[ipart+2*nparts]  = kp_;
    //     ( *delta )[ipart+0*nparts] = deltax;
    //     ( *delta )[ipart+1*nparts] = deltay;
    //     ( *delta )[ipart+2*nparts] = deltaz;
    // 
    // }
    
} // END Interpolator3D2Order


void InterpolatorAM2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc )
{
    // // Static cast of the electromagnetic fields
    // Field3D *Env_A_abs_3D = static_cast<Field3D *>( EMfields->Env_A_abs_ );
    // Field3D *Env_Chi_3D = static_cast<Field3D *>( EMfields->Env_Chi_ );
    // Field3D *Env_E_abs_3D = static_cast<Field3D *>( EMfields->Env_E_abs_ );
    // 
    // // Normalized particle position
    // double xpn = particles.position( 0, ipart )*dx_inv_;
    // double ypn = particles.position( 1, ipart )*dy_inv_;
    // double zpn = particles.position( 2, ipart )*dz_inv_;
    // 
    // 
    // // Indexes of the central nodes
    // ip_ = round( xpn );
    // jp_ = round( ypn );
    // kp_ = round( zpn );
    // 
    // 
    // // Declaration and calculation of the coefficient for interpolation
    // double delta2;
    // 
    // 
    // deltax   = xpn - ( double )ip_;
    // delta2  = deltax*deltax;
    // coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
    // coeffxp_[1] = 0.75 - delta2;
    // coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );
    // 
    // deltay   = ypn - ( double )jp_;
    // delta2  = deltay*deltay;
    // coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
    // coeffyp_[1] = 0.75 - delta2;
    // coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );
    // 
    // deltaz   = zpn - ( double )kp_;
    // delta2  = deltaz*deltaz;
    // coeffzp_[0] = 0.5 * ( delta2-deltaz+0.25 );
    // coeffzp_[1] = 0.75 - delta2;
    // coeffzp_[2] = 0.5 * ( delta2+deltaz+0.25 );
    // 
    // 
    // //!\todo CHECK if this is correct for both primal & dual grids !!!
    // // First index for summation
    // ip_ = ip_ - i_domain_begin;
    // jp_ = jp_ - j_domain_begin;
    // kp_ = kp_ - k_domain_begin;
    // 
    // // -------------------------
    // // Interpolation of Env_A_abs_^(p,p,p)
    // // -------------------------
    // *( Env_A_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_A_abs_3D, ip_, jp_, kp_ );
    // 
    // // -------------------------
    // // Interpolation of Env_Chi_^(p,p,p)
    // // -------------------------
    // *( Env_Chi_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_Chi_3D, ip_, jp_, kp_ );
    // 
    // // -------------------------
    // // Interpolation of Env_E_abs_^(p,p,p)
    // // -------------------------
    // *( Env_E_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_E_abs_3D, ip_, jp_, kp_ );
    
} // END Interpolator3D2Order
