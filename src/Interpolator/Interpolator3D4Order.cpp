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
Interpolator3D4Order::Interpolator3D4Order( Params &params, Patch *patch ) : Interpolator3D( params, patch )
{

    dx_inv_ = 1.0/params.cell_length[0];
    dy_inv_ = 1.0/params.cell_length[1];
    dz_inv_ = 1.0/params.cell_length[2];
    
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
    double xpn = particles.position( 0, ipart )*dx_inv_;
    double ypn = particles.position( 1, ipart )*dy_inv_;
    double zpn = particles.position( 2, ipart )*dz_inv_;
    // Calculate coeffs
    coeffs( xpn, ypn, zpn );
    
    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], &coeffzp_[2], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], &coeffzp_[2], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], &coeffzd_[2], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], &coeffzd_[2], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], &coeffzd_[2], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], &coeffzp_[2], Bz3D, id_, jd_, kp_ );
    
} // END Interpolator3D4Order


void Interpolator3D4Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    double xpn = particles.position( 0, ipart )*dx_inv_;
    double ypn = particles.position( 1, ipart )*dy_inv_;
    double zpn = particles.position( 2, ipart )*dz_inv_;
    // Calculate coeffs
    coeffs( xpn, ypn, zpn );
    
    int nparts( particles.size() );
    
    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], &coeffzp_[2], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], &coeffzp_[2], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp_[2], &coeffyp_[2], &coeffzd_[2], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp_[2], &coeffyd_[2], &coeffzd_[2], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffyp_[2], &coeffzd_[2], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], &coeffzp_[2], Bz3D, id_, jd_, kp_ );
    // Interpolation of Jx^(d,p,p)
    JLoc->x = compute( &coeffxd_[2], &coeffyp_[2], &coeffzp_[2], Jx3D, id_, jp_, kp_ );
    // Interpolation of Jy^(p,d,p)
    JLoc->y = compute( &coeffxp_[2], &coeffyd_[2], &coeffzp_[2], Jy3D, ip_, jd_, kp_ );
    // Interpolation of Jz^(p,p,d)
    JLoc->z = compute( &coeffxp_[2], &coeffyp_[2], &coeffzd_[2], Jz3D, ip_, jp_, kd_ );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = compute( &coeffxp_[2], &coeffyp_[2], &coeffzp_[2], Rho3D, ip_, jp_, kp_ );
}

// Interpolator on another field than the basic ones
void Interpolator3D4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
{
    Field3D *F = static_cast<Field3D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[2] : &coeffxp_[2];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[2] : &coeffyp_[2];
    double *coeffz = F->isDual( 2 ) ? &coeffzd_[2] : &coeffzp_[2];
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    int *j = F->isDual( 1 ) ? &jd_ : &jp_;
    int *k = F->isDual( 2 ) ? &kd_ : &kp_;
    
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*dx_inv_;
        double ypn = particles.position( 1, ipart )*dy_inv_;
        double zpn = particles.position( 2, ipart )*dz_inv_;
        coeffs( xpn, ypn, zpn );
        FieldLoc[ipart] = compute( coeffx, coeffy, coeffz, F, *i, *j, *k );
    }
}

void Interpolator3D4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    
    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );
        //Buffering of iol and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *iold )[ipart+2*nparts]  = kp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;
        ( *delta )[ipart+2*nparts] = deltaz;
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
    
        int npart_tot = particles.size();
        for( int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields( EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset );
        }
        
    }
}

void Interpolator3D4Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order


void Interpolator3D4Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order


void Interpolator3D4Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3D4Order
