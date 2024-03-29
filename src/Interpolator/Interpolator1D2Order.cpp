#include "Interpolator1D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "LaserEnvelope.h"


using namespace std;

Interpolator1D2Order::Interpolator1D2Order( Params &params, Patch *patch ) : Interpolator1D( patch )
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
    double xpn = particles.position( 0, ipart )*dx_inv_;
    // Calculate coeffs
    int idx_p[1], idx_d[1];
    double delta_p[1];
    double coeffxp[3];
    double coeffxd[3];
    coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );

    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffxd, Ex1D,   idx_d[0] );
    *( BLoc+1*nparts ) = compute( coeffxd, By1D_m, idx_d[0] );
    *( BLoc+2*nparts ) = compute( coeffxd, Bz1D_m, idx_d[0] );

    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffxp, Ey1D,   idx_p[0] );
    *( ELoc+2*nparts ) = compute( coeffxp, Ez1D,   idx_p[0] );
    *( BLoc+0*nparts ) = compute( coeffxp, Bx1D_m, idx_p[0] );

}//END Interpolator1D2Order

void Interpolator1D2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
{
    int ipart = *istart;

    double *ELoc = &( smpi->dynamics_Epart[ithread][ipart] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][ipart] );
    double *BLocyBTIS3;
    double *BLoczBTIS3;
    if(smpi->use_BTIS3){
        BLocyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][ipart] );
        BLoczBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][ipart] );
    }

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
    Field1D *By1DBTIS3;
    Field1D *Bz1DBTIS3;
    if (smpi->use_BTIS3){
        By1DBTIS3 = static_cast<Field1D *>( EMfields->By_mBTIS3 );
        Bz1DBTIS3 = static_cast<Field1D *>( EMfields->Bz_mBTIS3 );
    }

    // Particle position (in units of the spatial-step)
    double xpn = particles.position( 0, ipart )*dx_inv_;
    // Calculate coeffs
    int idx_p[1], idx_d[1];
    double delta_p[1];
    double coeffxp[3];
    double coeffxd[3];
    coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );

    int nparts( particles.numberOfParticles() );

    // Interpolate the fields from the Dual grid : Ex, By, Bz
    *( ELoc+0*nparts ) = compute( coeffxd, Ex1D,   idx_d[0] );
    *( BLoc+1*nparts ) = compute( coeffxd, By1D_m, idx_d[0] );
    *( BLoc+2*nparts ) = compute( coeffxd, Bz1D_m, idx_d[0] );

    // Interpolate the fields from the Primal grid : Ey, Ez, Bx
    *( ELoc+1*nparts ) = compute( coeffxp, Ey1D,   idx_p[0] );
    *( ELoc+2*nparts ) = compute( coeffxp, Ez1D,   idx_p[0] );
    *( BLoc+0*nparts ) = compute( coeffxp, Bx1D_m, idx_p[0] );

    // Interpolate the fields from the Primal grid : Jy, Jz, Rho
    JLoc->y = compute(     coeffxp, Jy1D,  idx_p[0] );
    JLoc->z = compute(     coeffxp, Jz1D,  idx_p[0] );
    ( *RhoLoc ) = compute( coeffxp, Rho1D, idx_p[0] );

    // Interpolate the fields from the Dual grid : Jx
    JLoc->x = compute( coeffxd, Jx1D,  idx_d[0] );
    
    if (smpi->use_BTIS3){
        *( BLocyBTIS3+0*nparts ) = compute( coeffxp, By1DBTIS3, idx_p[0] );
        *( BLoczBTIS3+0*nparts ) = compute( coeffxp, Bz1DBTIS3, idx_p[0] );
    }

}

// Interpolator on another field than the basic ones
void Interpolator1D2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field1D *F = static_cast<Field1D *>( *field );
    int idx_p[1], idx_d[1];
    double delta_p[1];
    double coeffxp[3];
    double coeffxd[3];
    double *coeff = F->isDual( 0 ) ? coeffxd : coeffxp;
    int *i = F->isDual( 0 ) ? &idx_d[0] : &idx_p[0];

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*dx_inv_;
        coeffs( xpn, idx_p, idx_d, coeffxp, coeffxd, delta_p );
        FieldLoc[ipart] = compute( coeff, F, *i );
    }
}

void Interpolator1D2Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
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
    int nparts = particles.numberOfParticles();
    
    if (!smpi->use_BTIS3){ // without BTIS-3 interpolation
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
            *( iold+0*nparts+ipart)   = idx_p[0];
            *( delta+0*nparts+ipart)  = delta_p[0];

        } // end ipart loop
    } else { // with B-TIS3 interpolation
      
      Field1D *By1D_mBTIS3 = static_cast<Field1D *>( EMfields->By_mBTIS3 );
      Field1D *Bz1D_mBTIS3 = static_cast<Field1D *>( EMfields->Bz_mBTIS3 );
      double  *BypartBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][0]  );
      double  *BzpartBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][0]  );
      
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
          // Interpolation of ByBTIS3^(p)
          *( BypartBTIS3+0*nparts )  = compute( coeffxp, By1D_mBTIS3, idx_p[0] );
          // Interpolation of BzBTIS3^(p)
          *( BzpartBTIS3+0*nparts )  = compute( coeffxp, Bz1D_mBTIS3, idx_p[0] );
          
          //Buffering of iol and delta
          *( iold+0*nparts+ipart)   = idx_p[0];
          *( delta+0*nparts+ipart)  = delta_p[0];

      } // end ipart loop
      
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

        int npart_tot = particles.numberOfParticles();
        for( int ipart=0 ; ipart<npart_tot; ipart++ ) {
            fields( EMfields, particles, ipart, offset, buffer+ipart, buffer+ipart+3*offset );
        }

    }
}

void Interpolator1D2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
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
    int nparts( particles.numberOfParticles() );

    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
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
            *( PHIpart+0*nparts+ipart ) = compute( coeffxp, Phi1D, idx_p[0] );
            // Interpolation of GradPhix^(p)
            *( GradPHIpart+0*nparts+ipart ) = compute( coeffxp, GradPhix1D, idx_p[0] );
            // Interpolation of GradPhiy^(p)
            *( GradPHIpart+1*nparts+ipart ) = compute( coeffxp, GradPhiy1D, idx_p[0] );
            // Interpolation of GradPhiz^(p)
            *( GradPHIpart+2*nparts+ipart ) = compute( coeffxp, GradPhiz1D, idx_p[0] );

            //Buffering of iol and delta
            *( iold+0*nparts+ipart)  = idx_p[0];
            *( delta+0*nparts+ipart) = delta_p[0];

        } // end ipart loop
    } else { // with B-TIS3 interpolation
        Field1D *By1D_mBTIS3 = static_cast<Field1D *>( EMfields->By_mBTIS3 );
        Field1D *Bz1D_mBTIS3 = static_cast<Field1D *>( EMfields->Bz_mBTIS3 );
        double  *BypartBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][0]  );
        double  *BzpartBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][0]  );
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
            // Interpolation of ByBTIS3^(p)
            *( BypartBTIS3+0*nparts )  = compute( coeffxp, By1D_mBTIS3, idx_p[0] );
            // Interpolation of BzBTIS3^(p)
            *( BzpartBTIS3+0*nparts )  = compute( coeffxp, Bz1D_mBTIS3, idx_p[0] );
            // Interpolation of Phi^(p)
            *( PHIpart+0*nparts+ipart )     = compute( coeffxp, Phi1D, idx_p[0] );
            // Interpolation of GradPhix^(p)
            *( GradPHIpart+0*nparts+ipart ) = compute( coeffxp, GradPhix1D, idx_p[0] );
            // Interpolation of GradPhiy^(p)
            *( GradPHIpart+1*nparts+ipart ) = compute( coeffxp, GradPhiy1D, idx_p[0] );
            // Interpolation of GradPhiz^(p)
            *( GradPHIpart+2*nparts+ipart ) = compute( coeffxp, GradPhiz1D, idx_p[0] );

            //Buffering of iol and delta
            *( iold+0*nparts+ipart)  = idx_p[0];
            *( delta+0*nparts+ipart) = delta_p[0];

        } // end ipart loop
    }

} // END Interpolator1D2Order

void Interpolator1D2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
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
    int nparts( particles.numberOfParticles() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;

        // Calculate coeffs

        int idx_p[1];
        double delta_p[1];
        double coeffxp[3];

        coeffs( xpn, idx_p, NULL, coeffxp, NULL, delta_p );

        // Interpolation of Phi^(p)
        *( PHI_mpart+0*nparts+ipart )     = compute( coeffxp, Phi_m1D, idx_p[0] );
        // Interpolation of GradPhix^(p)
        *( GradPHI_mpart+0*nparts+ipart ) = compute( coeffxp, GradPhix_m1D, idx_p[0] );
        // Interpolation of GradPhiy^(p)
        *( GradPHI_mpart+1*nparts+ipart ) = compute( coeffxp, GradPhiy_m1D, idx_p[0] );
        // Interpolation of GradPhiz^(p)
        *( GradPHI_mpart+2*nparts+ipart ) = compute( coeffxp, GradPhiz_m1D, idx_p[0] );

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
    int idx_p[1];
    double delta_p[1];
    double coeffxp[3];
    coeffs( xpn, idx_p, NULL, coeffxp, NULL, delta_p );

    // -------------------------
    // Interpolation of Env_A_abs_^(p)
    // -------------------------
    *( Env_A_abs_Loc )  = compute( coeffxp, Env_A_abs_1D, idx_p[0] ); 

    // -------------------------
    // Interpolation of Env_Chi_^(p)
    // -------------------------
    *( Env_Chi_Loc )    = compute( coeffxp, Env_Chi_1D, idx_p[0] ); 

    // -------------------------
    // Interpolation of Env_E_abs_^(p)
    // -------------------------
    *( Env_E_abs_Loc )  = compute( coeffxp, Env_E_abs_1D, idx_p[0] ); 

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( coeffxp, Env_Ex_abs_1D, idx_p[0] ); 

} // END Interpolator1D2Order

void Interpolator1D2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field1D *Env_Eabs = static_cast<Field1D *>( EMfields->Env_E_abs_ );

    std::vector<double> *Env_Eabs_part = &( smpi->dynamics_EnvEabs_part[ithread] );

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*dx_inv_;

        int idx_p[1];
        double delta_p[1];
        double coeffxp[3];
        coeffs( xpn, idx_p, NULL, coeffxp, NULL, delta_p );

        // ---------------------------------
        // Interpolation of Env_E_abs^(p)
        // ---------------------------------
        ( *Env_Eabs_part )[ipart] = compute( coeffxp, Env_Eabs, idx_p[0] );

        // In 1D the Env_Ex_abs field is always zero

    }

} // END Interpolator1D2Order
