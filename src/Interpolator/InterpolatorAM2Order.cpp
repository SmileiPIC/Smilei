#include "InterpolatorAM2Order.h"

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
