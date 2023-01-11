#include "Interpolator2DWT4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2DWT4Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2DWT4Order::Interpolator2DWT4Order( Params &params, Patch *patch ) : Interpolator2D( patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    dt_ov_D[0] = params.timestep/params.cell_length[0];
    dt_ov_D[1] = params.timestep/params.cell_length[1];
    dt2_ov_D2[0] = dt_ov_D[0]*dt_ov_D[0];
    dt2_ov_D2[1] = dt_ov_D[1]*dt_ov_D[1];
    D_ov_96dt[0] = 1.0/96.0/dt_ov_D[0];
    D_ov_96dt[1] = 1.0/96.0/dt_ov_D[1];


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
void Interpolator2DWT4Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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

    // Interpolation of Ex^(d,pt)
    *( ELoc+0*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(pt,d)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(pt,pt)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[2], &coeffypt_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(pt,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], By2D, id_, jp_ );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd_[2], &coeffyd_[2], Bz2D, id_, jd_ );
} // END Interpolator2DWT4Order

// -----------------------------------------------------------------------------
//
//! Interpolation of all fields and currents for a single particles
//! located at istart.
//! This version is not vectorized.
//! The input parameter iend not used for now, probes are interpolated one by one for now.
//
// -----------------------------------------------------------------------------
void Interpolator2DWT4Order::fieldsAndCurrents( ElectroMagn *EMfields,
                                                Particles &particles,
                                                SmileiMPI *smpi,
                                                int *istart,
                                                int *,
                                                int ithread,
                                                LocalFields *JLoc,
                                                double *RhoLoc )
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

    int nparts( particles.numberOfParticles() );

    // Interpolation of Ex^(d,pt)
    *( ELoc+0*nparts ) =  compute( &coeffxd_[2], &coeffypt_[2], Ex2D, id_, jp_ );
    // Interpolation of Ey^(pt,d)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], Ey2D, ip_, jd_ );
    // Interpolation of Ez^(pt,pt)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[2], &coeffypt_[2], Ez2D, ip_, jp_ );
    // Interpolation of Bx^(pt,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[2], &coeffyd_[2], Bx2D, ip_, jd_ );
    // Interpolation of By^(d,pt)
    *( BLoc+1*nparts ) = compute( &coeffxd_[2], &coeffypt_[2], By2D, id_, jp_ );
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
void Interpolator2DWT4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field2D *F = static_cast<Field2D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[2] : &coeffxpt_[2];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[2] : &coeffypt_[2];
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    int *j = F->isDual( 1 ) ? &jd_ : &jp_;

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        coeffs( xpn, ypn );
        FieldLoc[ipart] = compute( coeffx, coeffy, F, *i, *j );
    }
}

void Interpolator2DWT4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // std::cerr << "ipart: " << ipart
        //           << " x: " << particles.position( 0, ipart )
        //           << " y: " << particles.position( 1, ipart );

        //Interpolation on current particle
        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );
        //Buffering of iol and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;

        // std::cerr << " Ex: " << ( *Epart )[ipart+0*nparts]
        //           << " Ey: " << ( *Epart )[ipart+1*nparts]
        //           << " Ez: " << ( *Epart )[ipart+2*nparts]
        //           << " Bx: " << ( *Bpart )[ipart+0*nparts]
        //           << " By: " << ( *Bpart )[ipart+1*nparts]
        //           << " Bz: " << ( *Bpart )[ipart+2*nparts]
        //           << " iold: " << ( *iold )[ipart+0*nparts]
        //           << std::endl;

    }

}

// -----------------------------------------------------------------------------
//! Interpolator specific to tracked particles. A selection of particles may be provided
// -----------------------------------------------------------------------------
void Interpolator2DWT4Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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


void Interpolator2DWT4Order::fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator2DWT4Order::timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator2DWT4Order::envelopeAndSusceptibility( ElectroMagn *, Particles &, int, double *, double *, double *, double * )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
