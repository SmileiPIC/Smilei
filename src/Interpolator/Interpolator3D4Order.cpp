#include "Interpolator3D4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D4Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D4Order::Interpolator3D4Order( Params &params, Patch *patch ) : Interpolator3D( patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    d_inv_[2] = 1.0/params.cell_length[2];

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
// 4th Order Interpolation of the fields at a the particle position (5 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D4Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{
    // Static cast of the electromagnetic fields
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    double zpn = particles.position( 2, ipart )*d_inv_[2];
    // Compute coeffs
    int    idx_p[3], idx_d[3];
    double delta_p[3];
    double coeffxp[5], coeffyp[5], coeffzp[5];
    double coeffxd[5], coeffyd[5], coeffzd[5];
    coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[2], &coeffyp[2], &coeffzp[2], Ex3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp[2], &coeffyd[2], &coeffzp[2], Ey3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp[2], &coeffyp[2], &coeffzd[2], Ez3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[2], &coeffyd[2], &coeffzd[2], Bx3D, idx_p[0], idx_d[1], idx_d[2] );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd[2], &coeffyp[2], &coeffzd[2], By3D, idx_d[0], idx_p[1], idx_d[2] );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd[2], &coeffyd[2], &coeffzp[2], Bz3D, idx_d[0], idx_d[1], idx_p[2] );

} // END Interpolator3D4Order

void Interpolator3D4Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    // Interpolate E, B
    // Compute coefficient for ipart position
    // Static cast of the electromagnetic fields
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );
    Field3D *Jx3D = static_cast<Field3D *>( EMfields->Jx_ );
    Field3D *Jy3D = static_cast<Field3D *>( EMfields->Jy_ );
    Field3D *Jz3D = static_cast<Field3D *>( EMfields->Jz_ );
    Field3D *Rho3D= static_cast<Field3D *>( EMfields->rho_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    double zpn = particles.position( 2, ipart )*d_inv_[2];
    // Compute coeffs
    int    idx_p[3], idx_d[3];
    double delta_p[3];
    double coeffxp[5], coeffyp[5], coeffzp[5];
    double coeffxd[5], coeffyd[5], coeffzd[5];
    coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

    int nparts( particles.numberOfParticles() );

    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[2], &coeffyp[2], &coeffzp[2], Ex3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp[2], &coeffyd[2], &coeffzp[2], Ey3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp[2], &coeffyp[2], &coeffzd[2], Ez3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[2], &coeffyd[2], &coeffzd[2], Bx3D, idx_p[0], idx_d[1], idx_d[2] );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd[2], &coeffyp[2], &coeffzd[2], By3D, idx_d[0], idx_p[1], idx_d[2] );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd[2], &coeffyd[2], &coeffzp[2], Bz3D, idx_d[0], idx_d[1], idx_p[2] );
    // Interpolation of Jx^(d,p,p)
    JLoc->x = compute( &coeffxd[2], &coeffyp[2], &coeffzp[2], Jx3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Jy^(p,d,p)
    JLoc->y = compute( &coeffxp[2], &coeffyd[2], &coeffzp[2], Jy3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Jz^(p,p,d)
    JLoc->z = compute( &coeffxp[2], &coeffyp[2], &coeffzd[2], Jz3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = compute( &coeffxp[2], &coeffyp[2], &coeffzp[2], Rho3D, idx_p[0], idx_p[1], idx_p[2] );

}

// Interpolator on another field than the basic ones
void Interpolator3D4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field3D *F = static_cast<Field3D *>( *field );
    int    idx_p[3], idx_d[3];
    double delta_p[3];
    double coeffxp[5], coeffyp[5], coeffzp[5];
    double coeffxd[5], coeffyd[5], coeffzd[5];
    double *coeffx = F->isDual( 0 ) ? &coeffxd[2] : &coeffxp[2];
    double *coeffy = F->isDual( 1 ) ? &coeffyd[2] : &coeffyp[2];
    double *coeffz = F->isDual( 2 ) ? &coeffzd[2] : &coeffzp[2];
    int *i = F->isDual( 0 ) ? &idx_d[0] : &idx_p[0];
    int *j = F->isDual( 1 ) ? &idx_d[1] : &idx_p[1];
    int *k = F->isDual( 2 ) ? &idx_d[2] : &idx_p[2];

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];
        coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

        FieldLoc[ipart] = compute( coeffx, coeffy, coeffz, F, *i, *j, *k );
    }
}

void Interpolator3D4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
{
    double *const __restrict__ ELoc = &( smpi->dynamics_Epart[ithread][0] );
    double *const __restrict__ BLoc = &( smpi->dynamics_Bpart[ithread][0] );

    int    *const __restrict__ iold     = &( smpi->dynamics_iold[ithread][0] );
    double *const __restrict__ delta = &( smpi->dynamics_deltaold[ithread][0] );

    const double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    const double *const __restrict__ position_y = particles.getPtrPosition( 1 );
    const double *const __restrict__ position_z = particles.getPtrPosition( 2 );

    // Static cast of the electromagnetic fields
    Field3D *Ex3D = static_cast<Field3D *>( EMfields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( EMfields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( EMfields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( EMfields->Bx_m );
    Field3D *By3D = static_cast<Field3D *>( EMfields->By_m );
    Field3D *Bz3D = static_cast<Field3D *>( EMfields->Bz_m );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        // Coeffs
        int idx_p[3], idx_d[3];
        double delta_p[3];
        double coeffxp[5], coeffyp[5], coeffzp[5];
        double coeffxd[5], coeffyd[5], coeffzd[5];

        // Normalized particle position
        const double xpn = position_x[ ipart ]*d_inv_[0];
        const double ypn = position_y[ ipart ]*d_inv_[1];
        const double zpn = position_z[ ipart ]*d_inv_[2];

        coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

        // Interpolation of Ex^(d,p,p)
        *( ELoc+0*nparts+ipart ) = compute( &coeffxd[2], &coeffyp[2], &coeffzp[2], Ex3D, idx_d[0], idx_p[1], idx_p[2] );
        // Interpolation of Ey^(p,d,p)
        *( ELoc+1*nparts+ipart ) = compute( &coeffxp[2], &coeffyd[2], &coeffzp[2], Ey3D, idx_p[0], idx_d[1], idx_p[2] );
        // Interpolation of Ez^(p,p,d)
        *( ELoc+2*nparts+ipart ) = compute( &coeffxp[2], &coeffyp[2], &coeffzd[2], Ez3D, idx_p[0], idx_p[1], idx_d[2] );
        // Interpolation of Bx^(p,d,d)
        *( BLoc+0*nparts+ipart ) = compute( &coeffxp[2], &coeffyd[2], &coeffzd[2], Bx3D, idx_p[0], idx_d[1], idx_d[2] );
        // Interpolation of By^(d,p,d)
        *( BLoc+1*nparts+ipart ) = compute( &coeffxd[2], &coeffyp[2], &coeffzd[2], By3D, idx_d[0], idx_p[1], idx_d[2] );
        // Interpolation of Bz^(d,d,p)
        *( BLoc+2*nparts+ipart ) = compute( &coeffxd[2], &coeffyd[2], &coeffzp[2], Bz3D, idx_d[0], idx_d[1], idx_p[2] );

        //Buffering of iol and delta
        *( iold+0*nparts+ipart )  = idx_p[0];
        *( iold+1*nparts+ipart )  = idx_p[1];
        *( iold+2*nparts+ipart )  = idx_p[2];
        *( delta+0*nparts+ipart ) = delta_p[0];
        *( delta+1*nparts+ipart ) = delta_p[1];
        *( delta+2*nparts+ipart ) = delta_p[2];
    }
}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator3D4Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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

void Interpolator3D4Order::fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order


void Interpolator3D4Order::timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order


void Interpolator3D4Order::envelopeAndSusceptibility( ElectroMagn *, Particles &, int, double *, double *, double *, double * )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order
