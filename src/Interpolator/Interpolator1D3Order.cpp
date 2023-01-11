#include "Interpolator1D3Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "Patch.h"

using namespace std;

Interpolator1D3Order::Interpolator1D3Order( Params &params, Patch *patch ) : Interpolator1D( patch )
{
    dx_inv_   = 1.0/params.cell_length[0];

    dble_1ov6 = 1.0/6.0;
    dble_2ov3 = 2.0/3.0;
}

/***********************************************************************
    Interpolate the field fx defined on the primal grid
    with size nstp_x and space step stp_x_inv at the position
    xj and return the value fxj
***********************************************************************/
void Interpolator1D3Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
{
    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D *Ex1D     = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D     = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D     = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D_m   = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D_m   = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D_m   = static_cast<Field1D *>( EMfields->Bz_m );

    // Calculate the normalized positions
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

}//END Interpolator1D3Order


void Interpolator1D3Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
{

    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
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

    // Calculate the normalized positions
    double xjn = particles.position( 0, ipart )*dx_inv_;
    // Calculate coeffs
    coeffs( xjn );

    int nparts( particles.numberOfParticles() );

    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffd_, Ex1D,   id_ );
    *( BLoc+1*nparts ) = compute( coeffd_, By1D_m, id_ );
    *( BLoc+2*nparts ) = compute( coeffd_, Bz1D_m, id_ );

    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffp_, Ey1D,   ip_ );
    *( ELoc+2*nparts ) = compute( coeffp_, Ez1D,   ip_ );
    *( BLoc+0*nparts ) = compute( coeffp_, Bx1D_m, ip_ );

    // Primal Grid : Jy, Jz, Rho
    JLoc->y = compute( coeffp_, Jy1D,  ip_ );
    JLoc->z = compute( coeffp_, Jz1D,  ip_ );
    ( *RhoLoc ) = compute( coeffp_, Rho1D, ip_ );

    // Dual Grid : Jx
    JLoc->x = compute( coeffd_, Jx1D,  id_ );

}


// Interpolator on another field than the basic ones
void Interpolator1D3Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
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

void Interpolator1D3Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
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
    const int nparts = particles.numberOfParticles();

    for (int ipart=*istart; ipart < *iend; ipart++){

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;

        // Calculate coeffs
        int idx_p[1], idx_d[1];
        double delta_p[1];
        double coeffxp[4];
        double coeffxd[4];

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
void Interpolator1D3Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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
