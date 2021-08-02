#include "Interpolator1D2OrderV.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "LaserEnvelope.h"


using namespace std;

Interpolator1D2OrderV::Interpolator1D2OrderV( Params &params, Patch *patch ) : Interpolator1D( params, patch )
{
    dx_inv_ = 1.0/params.cell_length[0];
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2OrderV::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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

}//END Interpolator1D2OrderV

void Interpolator1D2OrderV::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, LocalFields *JLoc, double *RhoLoc )
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
void Interpolator1D2OrderV::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *l1, double *l2, double *l3 )
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

void Interpolator1D2OrderV::fieldsWrapper( ElectroMagn *EMfields, Particles &particles,
                                          SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    int *iold = &( smpi->dynamics_iold[ithread][0] );
    double *delta = &( smpi->dynamics_deltaold[ithread][0] );

    //int nparts( ( smpi->dynamics_invgf[ithread] ).size() );
    int nparts = particles.size();

    double * __restrict__ position_x = particles.getPtrPosition(0);

    double * __restrict__ Epart_x= &( smpi->dynamics_Epart[ithread][0*nparts] );
    double * __restrict__ Epart_y= &( smpi->dynamics_Epart[ithread][1*nparts] );
    double * __restrict__ Epart_z= &( smpi->dynamics_Epart[ithread][2*nparts] );

    double * __restrict__ Bpart_x= &( smpi->dynamics_Bpart[ithread][0*nparts] );
    double * __restrict__ Bpart_y= &( smpi->dynamics_Bpart[ithread][1*nparts] );
    double * __restrict__ Bpart_z= &( smpi->dynamics_Bpart[ithread][2*nparts] );

    // Static cast of the electromagnetic fields
    Field1D *Ex1D     = static_cast<Field1D *>( EMfields->Ex_ );
    Field1D *Ey1D     = static_cast<Field1D *>( EMfields->Ey_ );
    Field1D *Ez1D     = static_cast<Field1D *>( EMfields->Ez_ );
    Field1D *Bx1D_m   = static_cast<Field1D *>( EMfields->Bx_m );
    Field1D *By1D_m   = static_cast<Field1D *>( EMfields->By_m );
    Field1D *Bz1D_m   = static_cast<Field1D *>( EMfields->Bz_m );

    double * __restrict__ Ex = &Ex1D->data_[0]; //EMfields->Ex_->data();
    double * __restrict__ Ey = &Ey1D->data_[0]; // EMfields->Ey_->data();
    double * __restrict__ Ez = &Ez1D->data_[0]; //EMfields->Ez_->data();
    double * __restrict__ Bx = &Bx1D_m->data_[0]; //EMfields->Bx_m->data();
    double * __restrict__ By = &By1D_m->data_[0]; //EMfields->By_m->data();
    double * __restrict__ Bz = &Bz1D_m->data_[0]; //EMfields->Bz_m->data();

    int idx;  //dual index computed
    int ipx;  //prim index computed
    double xjn;
    double xjmxi2;

    double coeffd[3];
    double coeffp[3];

    #pragma omp simd private(xjn, xjmxi2, coeffd, coeffp, idx, ipx)
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        //Interpolation on current particle

        //double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
        //double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );

        // Particle position (in units of the spatial-step)
        xjn = position_x[ipart]*dx_inv_;
        //xjn = particles.position( 0, ipart )*dx_inv_;

        // Dual
        idx      = round( xjn+0.5 );      // index of the central point
        xjmxi  = xjn - ( double )idx +0.5; // normalized distance to the central node
        xjmxi2 = xjmxi*xjmxi;            // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffd[0] = 0.5 * ( xjmxi2-xjmxi+0.25 );
        coeffd[1] = ( 0.75-xjmxi2 );
        coeffd[2] = 0.5 * ( xjmxi2+xjmxi+0.25 );

        idx -= index_domain_begin;

        // Primal
        ipx      = round( xjn );    // index of the central point
        xjmxi  = xjn -( double )ipx; // normalized distance to the central node
        xjmxi2 = xjmxi*xjmxi;   // square of the normalized distance to the central node

        // 2nd order interpolation on 3 nodes
        coeffp[0] = 0.5 * ( xjmxi2-xjmxi+0.25 );
        coeffp[1] = ( 0.75-xjmxi2 );
        coeffp[2] = 0.5 * ( xjmxi2+xjmxi+0.25 );

        ipx -= index_domain_begin;

        // // Interpolate the fields from the Dual grid : Ex, By, Bz
        Epart_x[ipart] = coeffd[0] * Ex[idx-1]   + coeffd[1] * Ex[idx]   + coeffd[2] * Ex[idx+1];
        Bpart_y[ipart] = coeffd[0] * By[idx-1]   + coeffd[1] * By[idx]   + coeffd[2] * By[idx+1];
        Bpart_z[ipart] = coeffd[0] * Bz[idx-1]   + coeffd[1] * Bz[idx]   + coeffd[2] * Bz[idx+1];

        // Interpolate the fields from the Primal grid : Ey, Ez, Bx
        Epart_y[ipart] = coeffp[0] * Ey[ipx-1]   + coeffp[1] * Ey[ipx]   + coeffp[2] * Ey[ipx+1];
        Epart_z[ipart] = coeffp[0] * Ez[ipx-1]   + coeffp[1] * Ez[ipx]   + coeffp[2] * Ez[ipx+1];
        Bpart_x[ipart] = coeffp[0] * Bx[ipx-1]   + coeffp[1] * Bx[ipx]   + coeffp[2] * Bx[ipx+1];

        // Interpolate the fields from the Dual grid : Ex, By, Bz
        // *( ELoc+0*nparts ) = compute( coeffd_, Ex1D,   id_ );
        // *( BLoc+1*nparts ) = compute( coeffd_, By1D_m, id_ );
        // *( BLoc+2*nparts ) = compute( coeffd_, Bz1D_m, id_ );

        // Interpolate the fields from the Primal grid : Ey, Ez, Bx
        // *( ELoc+1*nparts ) = compute( coeffp_, Ey1D,   ip_ );
        // *( ELoc+2*nparts ) = compute( coeffp_, Ez1D,   ip_ );
        // *( BLoc+0*nparts ) = compute( coeffp_, Bx1D_m, ip_ );

        //Buffering of iol and delta
        iold[ipart] = ipx;
        delta[ipart] = xjmxi;
    }

}

// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator1D2OrderV::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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

void Interpolator1D2OrderV::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D *Phi1D = static_cast<Field1D *>( EMfields->envelope->Phi_ );
    Field1D *GradPhix1D = static_cast<Field1D *>( EMfields->envelope->GradPhix_ );
    Field1D *GradPhiy1D = static_cast<Field1D *>( EMfields->envelope->GradPhiy_ );
    Field1D *GradPhiz1D = static_cast<Field1D *>( EMfields->envelope->GradPhiz_ );

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );


        // -------------------------
        // Interpolation of Phi^(p,p)
        // -------------------------
        ( *PHIpart )[ipart] = compute( &coeffp_[1], Phi1D, ip_ );

        // -------------------------
        // Interpolation of GradPhix^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffp_[1], GradPhix1D, ip_ );

        // -------------------------
        // Interpolation of GradPhiy^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffp_[1], GradPhiy1D, ip_ );

        // -------------------------
        // Interpolation of GradPhiz^(p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffp_[1], GradPhiz1D, ip_ );


        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *delta )[ipart+0*nparts] = deltax;

    }


} // END Interpolator1D2OrderV


void Interpolator1D2OrderV::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D *Phi_m1D = static_cast<Field1D *>( EMfields->envelope->Phi_m );
    Field1D *GradPhix_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhix_m );
    Field1D *GradPhiy_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhiy_m );
    Field1D *GradPhiz_m1D = static_cast<Field1D *>( EMfields->envelope->GradPhiz_m );

    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.size() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

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
        // Interpolation of Phiold^(p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffp_[1], Phi_m1D, ip_ );

        // -------------------------
        // Interpolation of GradPhix_m^(p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffp_[1], GradPhix_m1D, ip_ );

        // -------------------------
        // Interpolation of GradPhiyold^(p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffp_[1], GradPhiy_m1D, ip_ );

        // -------------------------
        // Interpolation of GradPhizold^(p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffp_[1], GradPhiz_m1D, ip_ );

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *delta )[ipart+0*nparts] = deltax;


    }

} // END Interpolator1D2OrderV


void Interpolator1D2OrderV::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
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

} // END Interpolator1D2OrderV


void Interpolator1D2OrderV::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int ipart_ref )
{
    // Static cast of the envelope fields
    Field1D *Env_Eabs = static_cast<Field1D *>( EMfields->Env_E_abs_ );

    std::vector<double> *Env_Eabs_part = &( smpi->dynamics_EnvEabs_part[ithread] );

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

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

        // ---------------------------------
        // Interpolation of Env_E_abs^(p)
        // ---------------------------------
        ( *Env_Eabs_part )[ipart] = compute( &coeffp_[1], Env_Eabs, ip_ );

    }

} // END Interpolator1D2OrderV
