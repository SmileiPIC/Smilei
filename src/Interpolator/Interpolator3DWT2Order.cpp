#include "Interpolator3DWT2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "Particles.h"
#include "LaserEnvelope.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3DWT2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3DWT2Order::Interpolator3DWT2Order( Params &params, Patch *patch ) : Interpolator3D( patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    d_inv_[2] = 1.0/params.cell_length[2];
    dt_ov_D[0] = params.timestep/params.cell_length[0];
    dt_ov_D[1] = params.timestep/params.cell_length[1];
    dt_ov_D[2] = params.timestep/params.cell_length[2];
    D_ov_8dt[0] = params.cell_length[0]/params.timestep/8.0;
    D_ov_8dt[1] = params.cell_length[1]/params.timestep/8.0;
    D_ov_8dt[2] = params.cell_length[2]/params.timestep/8.0;
    
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order WT Interpolation of the fields at a the particle position
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3DWT2Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    *( ELoc+0*nparts ) = compute( &coeffxd_[1], &coeffypt_[1], &coeffzpt_[1], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(pt,d,pt)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[1], &coeffyd_[1], &coeffzpt_[1], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(pt,pt,d)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[1], &coeffypt_[1], &coeffzd_[1], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(pt,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[1], &coeffyd_[1], &coeffzd_[1], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,pt,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[1], &coeffypt_[1], &coeffzd_[1], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,pt)
    *( BLoc+2*nparts ) = compute( &coeffxd_[1], &coeffyd_[1], &coeffzpt_[1], Bz3D, id_, jd_, kp_ );
} // END Interpolator3DWT2Order

void Interpolator3DWT2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    *( ELoc+0*nparts ) = compute( &coeffxd_[1], &coeffypt_[1], &coeffzpt_[1], Ex3D, id_, jp_, kp_ );
    // Interpolation of Ey^(pt,d,pt)
    *( ELoc+1*nparts ) = compute( &coeffxpt_[1], &coeffyd_[1], &coeffzpt_[1], Ey3D, ip_, jd_, kp_ );
    // Interpolation of Ez^(pt,pt,d)
    *( ELoc+2*nparts ) = compute( &coeffxpt_[1], &coeffypt_[1], &coeffzd_[1], Ez3D, ip_, jp_, kd_ );
    // Interpolation of Bx^(pt,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxpt_[1], &coeffyd_[1], &coeffzd_[1], Bx3D, ip_, jd_, kd_ );
    // Interpolation of By^(d,pt,d)
    *( BLoc+1*nparts ) = compute( &coeffxd_[1], &coeffypt_[1], &coeffzd_[1], By3D, id_, jp_, kd_ );
    // Interpolation of Bz^(d,d,pt)
    *( BLoc+2*nparts ) = compute( &coeffxd_[1], &coeffyd_[1], &coeffzpt_[1], Bz3D, id_, jd_, kp_ );
    // Interpolation of Jx^(d,p,p)
    JLoc->x = compute( &coeffxd_[1], &coeffyp_[1], &coeffzp_[1], Jx3D, id_, jp_, kp_ );
    // Interpolation of Jy^(p,d,p)
    JLoc->y = compute( &coeffxp_[1], &coeffyd_[1], &coeffzp_[1], Jy3D, ip_, jd_, kp_ );
    // Interpolation of Jz^(p,p,d)
    JLoc->z = compute( &coeffxp_[1], &coeffyp_[1], &coeffzd_[1], Jz3D, ip_, jp_, kd_ );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Rho3D, ip_, jp_, kp_ );

}

// Interpolator on another field than the basic ones
void Interpolator3DWT2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field3D *F = static_cast<Field3D *>( *field );
    double *coeffx = F->isDual( 0 ) ? &coeffxd_[1] : &coeffxpt_[1];
    double *coeffy = F->isDual( 1 ) ? &coeffyd_[1] : &coeffypt_[1];
    double *coeffz = F->isDual( 2 ) ? &coeffzd_[1] : &coeffzpt_[1];
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

void Interpolator3DWT2Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
{
    double *ELoc = &( smpi->dynamics_Epart[ithread][0] );
    double *BLoc = &( smpi->dynamics_Bpart[ithread][0] );

    int*     iold = &( smpi->dynamics_iold    [ithread][0] );
    double* delta = &( smpi->dynamics_deltaold[ithread][0] );

    double* position_x = particles.getPtrPosition(0);
    double* position_y = particles.getPtrPosition(1);
    double* position_z = particles.getPtrPosition(2);

    // Static cast of the electromagnetic fields
    double* Ex3D = EMfields->Ex_->data_;
    double* Ey3D = EMfields->Ey_->data_;
    double* Ez3D = EMfields->Ez_->data_;
    double* Bx3D = EMfields->Bx_m->data_;
    double* By3D = EMfields->By_m->data_;
    double* Bz3D = EMfields->Bz_m->data_;

    int nx_p = EMfields->Bx_m->dims_[0];
    int ny_p = EMfields->By_m->dims_[1];
    int nz_p = EMfields->Bz_m->dims_[2];
    int nx_d = nx_p+1;
    int ny_d = ny_p+1;
    int nz_d = nz_p+1;

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        //Interpolation on current particle

        // Normalized particle position
        double xpn = position_x[ ipart ]*d_inv_[0];
        double ypn = position_y[ ipart ]*d_inv_[1];
        double zpn = position_z[ ipart ]*d_inv_[2];
        // Calculate coeffs

        int idx_p[3], idx_d[3];
        double delta_p[3];
        double coeffxpt[3], coeffypt[3], coeffzpt[3];
        double coeffxd[3], coeffyd[3], coeffzd[3];

        coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxpt, coeffypt, coeffzpt, coeffxd, coeffyd, coeffzd, delta_p );

        // Interpolation of Ex^(d,pt,pt)
        *( ELoc+0*nparts+ipart ) = compute( &coeffxd[1], &coeffypt[1], &coeffzpt[1], Ex3D, idx_d[0], idx_p[1], idx_p[2], nx_d, ny_p, nz_p );
        // Interpolation of Ey^(pt,d,pt)
        *( ELoc+1*nparts+ipart ) = compute( &coeffxpt[1], &coeffyd[1], &coeffzpt[1], Ey3D, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );
        // Interpolation of Ez^(pt,pt,d)
        *( ELoc+2*nparts+ipart ) = compute( &coeffxpt[1], &coeffypt[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );
        // Interpolation of Bx^(pt,d,d)
        *( BLoc+0*nparts+ipart ) = compute( &coeffxpt[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2], nx_p, ny_d, nz_d );
        // Interpolation of By^(d,pt,d)
        *( BLoc+1*nparts+ipart ) = compute( &coeffxd[1], &coeffypt[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2], nx_d, ny_p, nz_d );
        // Interpolation of Bz^(d,d,pt)
        *( BLoc+2*nparts+ipart ) = compute( &coeffxd[1], &coeffyd[1], &coeffzpt[1], Bz3D, idx_d[0], idx_d[1], idx_p[2], nx_d, ny_d, nz_p );

        //Buffering of iol and delta
        iold[ipart+0*nparts]  = idx_p[0];
        iold[ipart+1*nparts]  = idx_p[1];
        iold[ipart+2*nparts]  = idx_p[2];
        delta[ipart+0*nparts] = delta_p[0];
        delta[ipart+1*nparts] = delta_p[1];
        delta[ipart+2*nparts] = delta_p[2];
    }

}


// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator3DWT2Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, vector<unsigned int> *selection )
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


void Interpolator3DWT2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field3D *Phi3D = static_cast<Field3D *>( EMfields->envelope->Phi_ );
    Field3D *GradPhix3D = static_cast<Field3D *>( EMfields->envelope->GradPhix_ );
    Field3D *GradPhiy3D = static_cast<Field3D *>( EMfields->envelope->GradPhiy_ );
    Field3D *GradPhiz3D = static_cast<Field3D *>( EMfields->envelope->GradPhiz_ );

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        fields( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );


        // -------------------------
        // Interpolation of Phi^(p,p,p)
        // -------------------------
        ( *PHIpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Phi3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhix^(p,p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhix3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhiy^(p,p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiy3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhiz^(p,p,p)
        // -------------------------
        ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiz3D, ip_, jp_, kp_ );


        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *iold )[ipart+2*nparts]  = kp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;
        ( *delta )[ipart+2*nparts] = deltaz;

    }


} // END Interpolator3DWT2Order


void Interpolator3DWT2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field3D *Phi_m3D = static_cast<Field3D *>( EMfields->envelope->Phi_m );
    Field3D *GradPhix_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhix_m );
    Field3D *GradPhiy_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhiy_m );
    Field3D *GradPhiz_m3D = static_cast<Field3D *>( EMfields->envelope->GradPhiz_m );

    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];

        coeffs( xpn, ypn, zpn );

        // -------------------------
        // Interpolation of Phi_m^(p,p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Phi_m3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhix_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhix_m3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhiy_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiy_m3D, ip_, jp_, kp_ );

        // -------------------------
        // Interpolation of GradPhiz_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], GradPhiz_m3D, ip_, jp_, kp_ );

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = ip_;
        ( *iold )[ipart+1*nparts]  = jp_;
        ( *iold )[ipart+2*nparts]  = kp_;
        ( *delta )[ipart+0*nparts] = deltax;
        ( *delta )[ipart+1*nparts] = deltay;
        ( *delta )[ipart+2*nparts] = deltaz;

    }

} // END Interpolator3DWT2Order


void Interpolator3DWT2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    // Static cast of the electromagnetic fields
    Field3D *Env_A_abs_3D  = static_cast<Field3D *>( EMfields->Env_A_abs_ );
    Field3D *Env_Chi_3D    = static_cast<Field3D *>( EMfields->Env_Chi_ );
    Field3D *Env_E_abs_3D  = static_cast<Field3D *>( EMfields->Env_E_abs_ );
    Field3D *Env_Ex_abs_3D = static_cast<Field3D *>( EMfields->Env_Ex_abs_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    double zpn = particles.position( 2, ipart )*d_inv_[2];


    // Indexes of the central nodes
    ip_ = round( xpn );
    jp_ = round( ypn );
    kp_ = round( zpn );


    // Declaration and calculation of the coefficient for interpolation
    double delta2;


    deltax   = xpn - ( double )ip_;
    delta2  = deltax*deltax;
    coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
    coeffxp_[1] = 0.75 - delta2;
    coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );

    deltay   = ypn - ( double )jp_;
    delta2  = deltay*deltay;
    coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
    coeffyp_[1] = 0.75 - delta2;
    coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );

    deltaz   = zpn - ( double )kp_;
    delta2  = deltaz*deltaz;
    coeffzp_[0] = 0.5 * ( delta2-deltaz+0.25 );
    coeffzp_[1] = 0.75 - delta2;
    coeffzp_[2] = 0.5 * ( delta2+deltaz+0.25 );


    //!\todo CHECK if this is correct for both primal & dual grids !!!
    // First index for summation
    ip_ = ip_ - i_domain_begin;
    jp_ = jp_ - j_domain_begin;
    kp_ = kp_ - k_domain_begin;

    // -------------------------
    // Interpolation of Env_A_abs_^(p,p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_A_abs_3D, ip_, jp_, kp_ );

    // -------------------------
    // Interpolation of Env_Chi_^(p,p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_Chi_3D, ip_, jp_, kp_ );

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_E_abs_3D, ip_, jp_, kp_ );


    // -------------------------
    // Interpolation of Env_E_abs_^(p,p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], Env_Ex_abs_3D, ip_, jp_, kp_ );



} // END Interpolator3DWT2Order

void Interpolator3DWT2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field3D *EnvEabs = static_cast<Field3D *>( EMfields->Env_E_abs_ );
    Field3D *EnvExabs = static_cast<Field3D *>( EMfields->Env_Ex_abs_ );


    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );


    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];


        // Indexes of the central nodes
        ip_ = round( xpn );
        jp_ = round( ypn );
        kp_ = round( zpn );


        // Declaration and calculation of the coefficient for interpolation
        double delta2;


        deltax   = xpn - ( double )ip_;
        delta2  = deltax*deltax;
        coeffxp_[0] = 0.5 * ( delta2-deltax+0.25 );
        coeffxp_[1] = 0.75 - delta2;
        coeffxp_[2] = 0.5 * ( delta2+deltax+0.25 );

        deltay   = ypn - ( double )jp_;
        delta2  = deltay*deltay;
        coeffyp_[0] = 0.5 * ( delta2-deltay+0.25 );
        coeffyp_[1] = 0.75 - delta2;
        coeffyp_[2] = 0.5 * ( delta2+deltay+0.25 );

        deltaz   = zpn - ( double )kp_;
        delta2  = deltaz*deltaz;
        coeffzp_[0] = 0.5 * ( delta2-deltaz+0.25 );
        coeffzp_[1] = 0.75 - delta2;
        coeffzp_[2] = 0.5 * ( delta2+deltaz+0.25 );


        //!\todo CHECK if this is correct for both primal & dual grids !!!
        // First index for summation
        ip_ = ip_ - i_domain_begin;
        jp_ = jp_ - j_domain_begin;
        kp_ = kp_ - k_domain_begin;

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], EnvEabs, ip_, jp_, kp_ );

        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp_[1], &coeffyp_[1], &coeffzp_[1], EnvExabs, ip_, jp_, kp_ );

    }


} // END Interpolator3DWT2Order
