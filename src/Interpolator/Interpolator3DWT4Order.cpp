#include "Interpolator3DWT4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3DWT4Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3DWT4Order::Interpolator3DWT4Order( Params &params, Patch *patch ) : Interpolator3D( patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    d_inv_[2] = 1.0/params.cell_length[2];
    dt_ov_D[0] = params.timestep/params.cell_length[0]; 
    dt_ov_D[1] = params.timestep/params.cell_length[1]; 
    dt_ov_D[2] = params.timestep/params.cell_length[2]; 
    dt2_ov_D2[0] = dt_ov_D[0]*dt_ov_D[0];
    dt2_ov_D2[1] = dt_ov_D[1]*dt_ov_D[1];
    dt2_ov_D2[2] = dt_ov_D[2]*dt_ov_D[2];
    D_ov_96dt[0] = 1.0/96.0/dt_ov_D[0];
    D_ov_96dt[1] = 1.0/96.0/dt_ov_D[1];
    D_ov_96dt[2] = 1.0/96.0/dt_ov_D[2];

    //double defined for use in coefficients
    dble_1_ov_6 = 1.0/6.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_115_ov_192 = 115.0/192.0;
    
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order WT Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3DWT4Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    // Calculate coeffs
    coeffs( xpn, ypn, zpn );

    // Interpolation of Ex^(d,pt,pt)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], &coeffzpt_[2], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(pt,d,pt)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], &coeffzpt_[2], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(pt,pt,d)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[2], &coeffypt_[2], &coeffzd_[2], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(pt,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], &coeffzd_[2], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,pt,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], &coeffzd_[2], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,pt)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], &coeffzpt_[2], Bz3D, id_, jd_, kp_ );

} // END Interpolator3DWT4Order


void Interpolator3DWT4Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    // Calculate coeffs
    coeffs( xpn, ypn, zpn );

    int nparts( particles.numberOfParticles() );

    // Interpolation of Ex^(d,pt,pt)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], &coeffzpt_[2], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(pt,d,pt)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], &coeffzpt_[2], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(pt,pt,d)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[2], &coeffypt_[2], &coeffzd_[2], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(pt,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], &coeffzd_[2], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,pt,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], &coeffzd_[2], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,pt)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], &coeffzpt_[2], Bz3D, id_, jd_, kp_ );
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
void Interpolator3DWT4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field3D *F = static_cast<Field3D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[2] : &coeffxpt_[2];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[2] : &coeffypt_[2];
    double *coeffz = F->isDual( 2 ) ? &coeffzd_[2] : &coeffzpt_[2];
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    int *j = F->isDual( 1 ) ? &jd_ : &jp_;
    int *k = F->isDual( 2 ) ? &kd_ : &kp_;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];
        coeffs( xpn, ypn, zpn );
        FieldLoc[ipart] = compute( coeffx, coeffy, coeffz, F, *i, *j, *k );
    }
}

void Interpolator3DWT4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, 
                                            unsigned int, int )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
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
void Interpolator3DWT4Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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

void Interpolator3DWT4Order::fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3DWT4Order


void Interpolator3DWT4Order::timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3DWT4Order


void Interpolator3DWT4Order::envelopeAndSusceptibility( ElectroMagn *, Particles &, int , double *, double *, double *, double * )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
} // END Interpolator3DWT4Order
