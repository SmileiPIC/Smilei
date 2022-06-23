#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "LaserEnvelope.h"


using namespace std;

Interpolator1D2Order::Interpolator1D2Order( Params &params, Patch *patch ) : Interpolator1D( params, patch )
{
    dx_inv_ = 1.0/params.cell_length[0];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{
    // Static cast of the electromagnetic fields
    Field1D *Ex1D     = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D     = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D     = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D_m   = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D_m   = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D_m   = static_cast<Field1D *>( EMfields->Bz_m );

    // Particle position (in units of the spatial-step)
    double xjn = particles.position( 0, ipart )*dx_inv_;
    // Calculate coeffs
    coeffs( xjn );

    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffd_, Ex1D,   id_ );
    *( BLoc+1*nparts ) = compute( coeffd_, By1D_m, id_ );
    *( BLoc+2*nparts ) = compute( coeffd_, Bz1D_m, id_ );

    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffp_, Ey1D,   ip_ );
    *( ELoc+2*nparts ) = compute( coeffp_, Ez1D,   ip_ );
    *( BLoc+0*nparts ) = compute( coeffp_, Bx1D_m, ip_ );

}//END Interpolator1D2Order

void Interpolator1D2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    // Static cast of the electromagnetic fields
    Field1D *Ex1D     = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D     = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D     = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D_m   = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D_m   = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D_m   = static_cast<Field1D *>( EMfields->Bz_m );
    Field1D *Jx1D     = static_cast<Field1D *>( EMfields->Jx_ );
    Field1D *Jy1D     = static_cast<Field1D *>( EMfields->Jy_ );
    Field1D *Jz1D     = static_cast<Field1D *>( EMfields->Jz_ );
    Field1D *Rho1D    = static_cast<Field1D *>( EMfields->rho_ );

    // Particle position (in units of the spatial-step)
    double xjn = particles.position( 0, ipart )*dx_inv_;
    // Calculate coeffs
    coeffs( xjn );

    int nparts( particles.size() );

    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffd_, Ex1D,   id_ );
    *( BLoc+1*nparts ) = compute( coeffd_, By1D_m, id_ );
    *( BLoc+2*nparts ) = compute( coeffd_, Bz1D_m, id_ );

    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffp_, Ey1D,   ip_ );
    *( ELoc+2*nparts ) = compute( coeffp_, Ez1D,   ip_ );
    *( BLoc+0*nparts ) = compute( coeffp_, Bx1D_m, ip_ );

    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    JLoc->y = compute( coeffp_, Jy1D,  ip_ );
    JLoc->z = compute( coeffp_, Jz1D,  ip_ );
    ( *RhoLoc ) = compute( coeffp_, Rho1D, ip_ );

    // Interpolate the fields from the Dual grid : Jx
    JLoc->x = compute( coeffd_, Jx1D,  id_ );

}

// Interpolator on another field than the basic ones
void Interpolator1D2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    Field1D *F = static_cast<Field1D *>( *field );
    double *coeff = F->isDual( 0 ) ? coeffd_ : coeffp_;
    int *i = F->isDual( 0 ) ? &id_ : &ip_;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xjn = particles.position( 0, ipart )*dx_inv_;
        coeffs( xjn );
        FieldLoc[ipart] = compute( coeff, F, *i );
    }
}

void Interpolator1D2Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell, int ipart_ref )
{
    double *Epart = &( smpi->dynamics_Epart[ithread][0] );
    double *Bpart = &( smpi->dynamics_Bpart[ithread][0] );
    int    *iold  = &( smpi->dynamics_iold[ithread][0] );
    double *delta = &( smpi->dynamics_deltaold[ithread][0] );

    // Static cast of the electromagnetic fields
    Field1D *Ex1D = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D = static_cast<Field1D *>( EMfields->Bz_m );


    //Loop on bin particles
    int nparts = particles.size();

    for (int ipart=*istart; ipart < *iend; ipart++){

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;
       
        // Calculate coeffs
        int idx_p[1], idx_d[1];
        double delta_p[1];
        double coeffxp[3];
        double coeffxd[3];

        coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );

        // Interpolation of Ex^(d)
        *( Epart+0*nparts+ipart ) = compute( coeffxd, Ex1D, idx_d[0] );
        // Interpolation of Ey^(p)
        *( Epart+1*nparts+ipart ) = compute( coeffxp, Ey1D, idx_p[0] );
        // Interpolation of Ez^(p)
        *( Epart+2*nparts+ipart ) = compute( coeffxp, Ez1D, idx_p[0] );
        // Interpolation of Bx^(p)
        *( Bpart+0*nparts+ipart ) = compute( coeffxp, Bx1D, idx_p[0] );
        // Interpolation of By^(d)
        *( Bpart+1*nparts+ipart ) = compute( coeffxd, By1D, idx_d[0] );
        // Interpolation of Bz^(d)
        *( Bpart+2*nparts+ipart ) = compute( coeffxd, Bz1D, idx_d[0] );

        //Buffering of iol and delta
        *( iold+0*nparts+ipart)  = idx_p[0];
        *( delta+0*nparts+ipart) = delta_p[0];

    }
}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator1D2Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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

void Interpolator1D2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    
    double *Epart       = &( smpi->dynamics_Epart[ithread][0] );
    double *Bpart       = &( smpi->dynamics_Bpart[ithread][0] );
    double *PHIpart     = &( smpi->dynamics_PHIpart[ithread][0] );
    double *GradPHIpart = &( smpi->dynamics_GradPHIpart[ithread][0] );
    int    *iold        = &( smpi->dynamics_iold[ithread][0] );
    double *delta       = &( smpi->dynamics_deltaold[ithread][0] );

    // Static cast of the electromagnetic fields
    Field1D *Ex1D = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D = static_cast<Field1D *>( EMfields->Bz_m );
    Field1D *Phi1D = static_cast<Field1D *>( EMfields->envelope->Phi_ );
    Field1D *GradPhix1D = static_cast<Field1D *>( EMfields->envelope->GradPhix_ );
    Field1D *GradPhiy1D = static_cast<Field1D *>( EMfields->envelope->GradPhiy_ );
    Field1D *GradPhiz1D = static_cast<Field1D *>( EMfields->envelope->GradPhiz_ );


    //Loop on bin particles
    int nparts = particles.size();

    for (int ipart=*istart; ipart < *iend; ipart++){

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;
       
        // Calculate coeffs
        int idx_p[1], idx_d[1];
        double delta_p[1];
        double coeffxp[3];
        double coeffxd[3];

        coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );

        // Interpolation of Ex^(d)
        *( Epart+0*nparts+ipart ) = compute( coeffxd, Ex1D, idx_d[0] );
        // Interpolation of Ey^(p)
        *( Epart+1*nparts+ipart ) = compute( coeffxp, Ey1D, idx_p[0] );
        // Interpolation of Ez^(p)
        *( Epart+2*nparts+ipart ) = compute( coeffxp, Ez1D, idx_p[0] );
        // Interpolation of Bx^(p)
        *( Bpart+0*nparts+ipart ) = compute( coeffxp, Bx1D, idx_p[0] );
        // Interpolation of By^(d)
        *( Bpart+1*nparts+ipart ) = compute( coeffxd, By1D, idx_d[0] );
        // Interpolation of Bz^(d)
        *( Bpart+2*nparts+ipart ) = compute( coeffxd, Bz1D, idx_d[0] );
        // Interpolation of Phi^(p)
        *( PHIpart+0*nparts+ipart ) = compute( coeffxp, Phi1D, idx_d[0] );
        // Interpolation of GradPhix^(p)
        *( GradPHIpart+0*nparts+ipart ) = compute( coeffxp, GradPhix1D, idx_d[0] );
        // Interpolation of GradPhiy^(p)
        *( GradPHIpart+1*nparts+ipart ) = compute( coeffxp, GradPhiy1D, idx_d[0] );
        // Interpolation of GradPhiz^(p)
        *( GradPHIpart+2*nparts+ipart ) = compute( coeffxp, GradPhiz1D, idx_d[0] );

        //Buffering of iol and delta
        *( iold+0*nparts+ipart)  = idx_p[0];
        *( delta+0*nparts+ipart) = delta_p[0];

    }

} // END Interpolator1D2Order

void Interpolator1D2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D *Phi_m1D = static_cast<Field1D *>( EMfields->envelope->Phi_m );
    Field1D *GradPhix_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhix_m );
    Field1D *GradPhiy_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhiy_m );
    Field1D *GradPhiz_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhiz_m );
    
    double *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread][0] );
    double *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread][0] );
    
    int    *iold  = &( smpi->dynamics_iold[ithread][0] );
    double *delta = &( smpi->dynamics_deltaold[ithread][0] );
    
    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
    
        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;
       
        // Calculate coeffs

        int idx_p[1], idx_d[1];
        double delta_p[1];
        double coeffxp[3];
        double coeffxd[3];

        coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );
        
        // Interpolation of Phi^(p)
        *( PHI_mpart+0*nparts+ipart ) = compute( coeffxp, Phi_m1D, idx_d[0] );
        // Interpolation of GradPhix^(p)
        *( GradPHI_mpart+0*nparts+ipart ) = compute( coeffxp, GradPhix_m1D, idx_d[0] );
        // Interpolation of GradPhiy^(p)
        *( GradPHI_mpart+1*nparts+ipart ) = compute( coeffxp, GradPhiy_m1D, idx_d[0] );
        // Interpolation of GradPhiz^(p)
        *( GradPHI_mpart+2*nparts+ipart ) = compute( coeffxp, GradPhiz_m1D, idx_d[0] );
        
        //Buffering of iol and delta
        *( iold+ipart+0*nparts)  = idx_p[0];
        *( delta+ipart+0*nparts) = delta_p[0];
        
        
    }
    
} // END Interpolator1D2Order


void Interpolator1D2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    // Static cast of the electromagnetic fields
    Field1D *Env_A_abs_1D  = static_cast<Field1D *>( EMfields->Env_A_abs_ );
    Field1D *Env_Chi_1D    = static_cast<Field1D *>( EMfields->Env_Chi_ );
    Field1D *Env_E_abs_1D  = static_cast<Field1D *>( EMfields->Env_E_abs_ );
    Field1D *Env_Ex_abs_1D = static_cast<Field1D *>( EMfields->Env_Ex_abs_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*dx_inv_;

    // Indexes of the central nodes
    ip_ = round( xpn );

    // Declaration and calculation of the coefficient for interpolation
    double deltax, delta2;


    deltax   = xpn - ( double )ip_;
    delta2  = deltax*deltax;
    coeffp_[0] = 0.5 * ( delta2-deltax+0.25 );
    coeffp_[1] = 0.75 - delta2;
    coeffp_[2] = 0.5 * ( delta2+deltax+0.25 );


    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - index_domain_begin;

    // -------------------------
    // Interpolation of Env_A_abs_^(p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffp_[1], Env_A_abs_1D, ip_ );

    // -------------------------
    // Interpolation of Env_Chi_^(p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffp_[1], Env_Chi_1D, ip_ );

    // -------------------------
    // Interpolation of Env_E_abs_^(p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffp_[1], Env_E_abs_1D, ip_ );

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffp_[1], Env_Ex_abs_1D, ip_ );

} // END Interpolator1D2Order

void Interpolator1D2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D *Env_Eabs = static_cast<Field1D *>( EMfields->Env_E_abs_ );
    
    std::vector<double> *Env_Eabs_part = &( smpi->dynamics_EnvEabs_part[ithread] );
    
    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        int idx_p[1];
        double delta_p[1];
        double coeffxp[3];

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;

        double delta2;
        
        // Primal
        idx_p[0]    = round( xpn );                 // index of the central point
        delta_p[0]  = xpn -( double )idx_p[0];      // normalized distance to the central node
        delta2      = pow( delta_p[0], 2 );         // square of the normalized distance to the central node
        
        // 2nd order interpolation on 3 nodes
        coeffxp[0]   = 0.5 * ( delta2-delta_p[0]+0.25 );
        coeffxp[1]   = ( 0.75-delta2 );
        coeffxp[2]   = 0.5 * ( delta2+delta_p[0]+0.25 );
        
        idx_p[0]   -= index_domain_begin;
    
        // ---------------------------------
        // Interpolation of Env_E_abs^(p)
        // ---------------------------------
        ( *Env_Eabs_part )[ipart] = compute( coeffxp, Env_Eabs, idx_p[0] );

        // In 1D the Env_Ex_abs field is always zero 
        
    }    
    
} // END Interpolator1D2Order
