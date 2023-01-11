#include "Interpolator1DWT4Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"

using namespace std;

Interpolator1DWT4Order::Interpolator1DWT4Order( Params &params, Patch *patch ) : Interpolator1D( patch )
{
    dx_inv_ = 1.0/params.cell_length[0];
    dt_ov_dx = params.timestep/params.cell_length[0]; 
    dt2_ov_dx2 = dt_ov_dx*dt_ov_dx;
    dx_ov_96dt = 1.0/96.0/dt_ov_dx;
    
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
void Interpolator1DWT4Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    
    // Interpolate the fields from the Primal grid (WT) : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffpt_, Ey1D,   ip_ );
    *( ELoc+2*nparts ) = compute( coeffpt_, Ez1D,   ip_ );
    *( BLoc+0*nparts ) = compute( coeffpt_, Bx1D_m, ip_ );
    
}//END Interpolator1DWT4Order

void Interpolator1DWT4Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    
    int nparts( particles.numberOfParticles() );
    
    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffd_, Ex1D,   id_ );
    *( BLoc+1*nparts ) = compute( coeffd_, By1D_m, id_ );
    *( BLoc+2*nparts ) = compute( coeffd_, Bz1D_m, id_ );
    
    // Interpolate the fields from the Primal (WT) grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffpt_, Ey1D,   ip_ );
    *( ELoc+2*nparts ) = compute( coeffpt_, Ez1D,   ip_ );
    *( BLoc+0*nparts ) = compute( coeffpt_, Bx1D_m, ip_ );
    
    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    JLoc->y = compute( coeffp_, Jy1D,  ip_ );
    JLoc->z = compute( coeffp_, Jz1D,  ip_ );
    ( *RhoLoc ) = compute( coeffp_, Rho1D, ip_ );
    
    // Interpolate the fields from the Dual grid : Jx
    JLoc->x = compute( coeffd_, Jx1D,  id_ );
    
}

// Interpolator on another field than the basic ones
void Interpolator1DWT4Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field1D *F = static_cast<Field1D *>( *field );
    double *coeff = F->isDual( 0 ) ? coeffd_ : coeffpt_;
    int *i = F->isDual( 0 ) ? &id_ : &ip_;
    
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xjn = particles.position( 0, ipart )*dx_inv_;
        coeffs( xjn );
        FieldLoc[ipart] = compute( coeff, F, *i );
    }
}

void Interpolator1DWT4Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<int> *iold = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );
    
    //Loop on bin particles
    int npart_tot = particles.numberOfParticles();
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle
        fields( EMfields, particles, ipart, npart_tot, &( *Epart )[ipart], &( *Bpart )[ipart] );
        //Buffering of iol and delta
        ( *iold )[ipart] = ip_;
        ( *delta )[ipart] = xjmxi;
    }
    
}


// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator1DWT4Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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

void Interpolator1DWT4Order::fieldsAndEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}


void Interpolator1DWT4Order::timeCenteredEnvelope( ElectroMagn *, Particles &, SmileiMPI *, int *, int *, int, int )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}

// probes like diagnostic !
void Interpolator1DWT4Order::envelopeAndSusceptibility( ElectroMagn *, Particles &, int , double *, double *, double *, double * )
{
    ERROR( "Projection and interpolation for the envelope model are implemented only for interpolation_order = 2" );
}
