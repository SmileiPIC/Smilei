#include "Interpolator2D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D4Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D4Order::Interpolator2D4Order( Params &params, Patch *patch ) : Interpolator2D( params, patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];


    //double defined for use in coefficients
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_48 = 1.0/48.0;
    dble_1_ov_16 = 1.0/16.0;
    dble_1_ov_12 = 1.0/12.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_1_ov_4 = 1.0/4.0;
    dble_1_ov_6 = 1.0/6.0;
    dble_115_ov_192 = 115.0/192.0;
    dble_5_ov_8 = 5.0/8.0;
    
}


// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D4Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    // Calculate coeffs
    coeffs( xpn, ypn );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], Bz2D, id_, jd_ );
} // END Interpolator2D4Order

void Interpolator2D4Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
{

    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );
    Field2D *Jx2D = static_cast<Field2D *>( EMfields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( EMfields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( EMfields->Jz_ );
    Field2D *Rho2D= static_cast<Field2D *>( EMfields->rho_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    // Calculate coeffs
    coeffs( xpn, ypn );

    int nparts( particles.size() );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) =  compute( &coeffxd_[2], &coeffyp_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], Bz2D, id_, jd_ );
    // Interpolation of Jx^(d,p)
    JLoc->x = compute( &coeffxd_[2], &coeffyp_[2], Jx2D, id_, jp_ );
    // Interpolation of Ey^(p,d)
    JLoc->y = compute( &coeffxp_[2], &coeffyd_[2], Jy2D, ip_, jd_ );
    // Interpolation of Ez^(p,p)
    JLoc->z = compute( &coeffxp_[2], &coeffyp_[2], Jz2D, ip_, jp_ );
    // Interpolation of Rho^(p,p)
    ( *RhoLoc ) = compute( &coeffxp_[2], &coeffyp_[2], Rho2D, ip_, jp_ );
}

//! Interpolator on another field than the basic ones
void Interpolator2D4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    Field2D *F = static_cast<Field2D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[2] : &coeffxp_[2];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[2] : &coeffyp_[2];
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    int *j = F->isDual( 1 ) ? &jd_ : &jp_;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        coeffs( xpn, ypn );
        FieldLoc[ipart] = compute( coeffx, coeffy, F, *i, *j );
    }
}

void Interpolator2D4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int scell, int ipart_ref )
{
    double *Epart = &( smpi->dynamics_Epart[ithread][0] );
    double *Bpart = &( smpi->dynamics_Bpart[ithread][0] );
    int *iold = &( smpi->dynamics_iold[ithread][0] );
    double *delta = &( smpi->dynamics_deltaold[ithread][0] );

    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    //Loop on bin particles
    int nparts( particles.size() );

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];

        // Coeffs
        int idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[5], coeffyp[5];
        double coeffxd[5], coeffyd[5];

        coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

        // Interpolation of Ex^(d,p)
        *( Epart+0*nparts ) = compute( &coeffxd[2], &coeffyp[2], Ex2D, idx_d[0], idx_p[1] );
        // Interpolation of Ey^(p,d)
        *( Epart+1*nparts ) = compute( &coeffxp[2], &coeffyd[2], Ey2D, idx_p[0], idx_d[1] );
        // Interpolation of Ez^(p,p)
        *( Epart+2*nparts ) = compute( &coeffxp[2], &coeffyp[2], Ez2D, idx_p[0], idx_p[1] );
        // Interpolation of Bx^(p,d)
        *( Bpart+0*nparts ) = compute( &coeffxp[2], &coeffyd[2], Bx2D, idx_p[0], idx_d[1] );
        // Interpolation of By^(d,p)
        *( Bpart+1*nparts ) = compute( &coeffxd[2], &coeffyp[2], By2D, idx_d[0], idx_p[1] );
        // Interpolation of Bz^(d,d)
        *( Bpart+2*nparts ) = compute( &coeffxd[2], &coeffyd[2], Bz2D, idx_d[0], idx_d[1] );

        //Buffering of iol and delta
        *( iold+0*nparts)  = idx_p[0];
        *( iold+1*nparts)  = idx_p[1];
        *( delta+0*nparts) = delta_p[0];
        *( delta+1*nparts) = delta_p[1];

    }

    
}

// -----------------------------------------------------------------------------
//! Interpolator specific to tracked particles. A selection of particles may be provided
// -----------------------------------------------------------------------------
void Interpolator2D4Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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


void Interpolator2D4Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator2D4Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator2D4Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
