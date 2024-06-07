#include "Interpolator2D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particles.h"
#include "LaserEnvelope.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator2D2Order::Interpolator2D2Order( Params &params, Patch *patch ) : Interpolator2D( patch )
{
    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    int idx_p[2], idx_d[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[3];
    double coeffxd[3], coeffyd[3];
    coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1] );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1] );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1] );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1] );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1] );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1] );
} // END Interpolator2D2Order

void Interpolator2D2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    Field2D *By2DBTIS3;
    Field2D *Bz2DBTIS3;
    if (smpi->use_BTIS3){
        By2DBTIS3 = static_cast<Field2D *>( EMfields->By_mBTIS3 );
        Bz2DBTIS3 = static_cast<Field2D *>( EMfields->Bz_mBTIS3 );
    }

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    // Calculate coeffs
    int idx_p[2], idx_d[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[3];
    double coeffxd[3], coeffyd[3];
    coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

    int nparts( particles.numberOfParticles() );

    // Interpolation of Ex^(d,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1] );
    // Interpolation of Ey^(p,d)
    *( ELoc+1*nparts ) = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1] );
    // Interpolation of Ez^(p,p)
    *( ELoc+2*nparts ) = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1] );
    // Interpolation of Bx^(p,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1] );
    // Interpolation of By^(d,p)
    *( BLoc+1*nparts ) = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1] );
    // Interpolation of Bz^(d,d)
    *( BLoc+2*nparts ) = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1] );
    // Interpolation of Jx^(d,p)
    JLoc->x = compute( &coeffxd[1], &coeffyp[1], Jx2D, idx_d[0], idx_p[1] );
    // Interpolation of Jy^(p,d)
    JLoc->y = compute( &coeffxp[1], &coeffyd[1], Jy2D, idx_p[0], idx_d[1] );
    // Interpolation of Jz^(p,p)
    JLoc->z = compute( &coeffxp[1], &coeffyp[1], Jz2D, idx_p[0], idx_p[1] );
    // Interpolation of Rho^(p,p)
    ( *RhoLoc ) = compute( &coeffxp[1], &coeffyp[1], Rho2D, idx_p[0], idx_p[1] );
    
    if (smpi->use_BTIS3){
        // Interpolation of ByBTIS3^(p,p)
        *( BLocyBTIS3+0*nparts ) = compute( &coeffxp[1], &coeffyp[1], By2DBTIS3, idx_p[0], idx_p[1] );
        // Interpolation of BzBTIS3^(p,d)
        *( BLoczBTIS3+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], Bz2DBTIS3, idx_p[0], idx_d[1] );
    }
}

//! Interpolator on another field than the basic ones
void Interpolator2D2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field2D *F = static_cast<Field2D *>( *field );
    int idx_p[2], idx_d[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[3];
    double coeffxd[3], coeffyd[3];
    double *coeffx = F->isDual( 0 ) ? &coeffxd[1] : &coeffxp[1];
    double *coeffy = F->isDual( 1 ) ? &coeffyd[1] : &coeffyp[1];
    int *i = F->isDual( 0 ) ? &idx_d[0] : &idx_p[0];
    int *j = F->isDual( 1 ) ? &idx_d[1] : &idx_p[1];

    for( int ipart=*istart ; ipart<*iend; ipart++ ) {
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );
        FieldLoc[ipart] = compute( coeffx, coeffy, F, *i, *j );
    }
}

// -----------------------------------------------------------------------------
//! Wrapper called by the particle dynamics section
// -----------------------------------------------------------------------------
void Interpolator2D2Order::fieldsWrapper(   ElectroMagn *EMfields,
                                            Particles &particles,
                                            SmileiMPI *smpi,
                                            int *istart,
                                            int *iend,
                                            int ithread,
                                            unsigned int,
                                            int )
{
    double *const __restrict__ ELoc  = smpi->dynamics_Epart[ithread].data();
    double *const __restrict__ BLoc  = smpi->dynamics_Bpart[ithread].data();

    int    *const __restrict__ iold  = smpi->dynamics_iold[ithread].data();
    double *const __restrict__ delta = smpi->dynamics_deltaold[ithread].data();

    const double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    const double *const __restrict__ position_y = particles.getPtrPosition( 1 );

    const double *const __restrict__ Ex2D = static_cast<Field2D *>( EMfields->Ex_ )->data();
    const double *const __restrict__ Ey2D = static_cast<Field2D *>( EMfields->Ey_ )->data();
    const double *const __restrict__ Ez2D = static_cast<Field2D *>( EMfields->Ez_ )->data();
    const double *const __restrict__ Bx2D = static_cast<Field2D *>( EMfields->Bx_m )->data();
    const double *const __restrict__ By2D = static_cast<Field2D *>( EMfields->By_m )->data();
    const double *const __restrict__ Bz2D = static_cast<Field2D *>( EMfields->Bz_m )->data();

#if defined(SMILEI_ACCELERATOR_GPU_OACC)    
    const int sizeofEx = EMfields->Ex_->size();
    const int sizeofEy = EMfields->Ey_->size();
    const int sizeofEz = EMfields->Ez_->size();
    const int sizeofBx = EMfields->Bx_m->size();
    const int sizeofBy = EMfields->By_m->size();
    const int sizeofBz = EMfields->Bz_m->size();
#endif

    const int ny_p = EMfields->By_m->dims_[1]; // primary_grid_size_in_y
    const int ny_d = ny_p + 1;                 // dual_grid_size_in_y

    // Loop on bin particles
    const int nparts = particles.numberOfParticles();

    const int first_index = *istart;
    const int last_index  = *iend;

    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
#if defined( SMILEI_ACCELERATOR_GPU_OMP )

    #pragma omp target map( to                                                     \
                            : i_domain_begin, j_domain_begin )                     \
        is_device_ptr /* map */ ( /* to: */                                        \
                                  position_x /* [first_index:npart_range_size] */, \
                                  position_y /* [first_index:npart_range_size] */ )
    #pragma omp teams distribute parallel for
#elif defined(SMILEI_ACCELERATOR_GPU_OACC)
    #pragma acc enter data create(this)
    #pragma acc update device(this)
    size_t interpolation_range_size = ( last_index + 1 * nparts ) - first_index;
    #pragma acc parallel present(ELoc [first_index:interpolation_range_size],\
                                 BLoc [first_index:interpolation_range_size],\
                                 iold [first_index:interpolation_range_size],\
                                 delta [first_index:interpolation_range_size],\
                                 Ex2D [0:sizeofEx],\
                                 Ey2D [0:sizeofEy],\
                                 Ez2D [0:sizeofEz],\
                                 Bx2D [0:sizeofBx],\
                                 By2D [0:sizeofBy],\
                                 Bz2D [0:sizeofBz])\
    deviceptr(position_x, position_y)              \
    copyin(d_inv_[0:2])
    #pragma acc loop gang worker vector
#endif
    for( int ipart = first_index; ipart < last_index; ipart++ ) {
        // Interpolation on current particle

        // Normalized particle position
        const double xpn = position_x[ipart] * d_inv_[0];
        const double ypn = position_y[ipart] * d_inv_[1];

        // Calculate coeffs
        int    idx_p[2], idx_d[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        double coeffxd[3], coeffyd[3];

        coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

        // Interpolation of Ex^(d,p)
        ELoc[0*nparts+ipart] = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1], ny_p );
        // Interpolation of Ey^(p,d)
        ELoc[1*nparts+ipart] = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1], ny_d );
        // Interpolation of Ez^(p,p)
        ELoc[2*nparts+ipart] = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1], ny_p );
        // Interpolation of Bx^(p,d)
        BLoc[0*nparts+ipart] = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1], ny_d );
        // Interpolation of By^(d,p)
        BLoc[1*nparts+ipart] = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1], ny_p );
        // Interpolation of Bz^(d,d)
        BLoc[2*nparts+ipart] = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1], ny_d );

        // Buffering of iol and delta
        iold[0*nparts+ipart]  = idx_p[0];
        iold[1*nparts+ipart]  = idx_p[1];
        delta[0*nparts+ipart] = delta_p[0];
        delta[1*nparts+ipart] = delta_p[1];
        
    }
    #if defined(SMILEI_ACCELERATOR_GPU_OACC)
        #pragma acc exit data delete(this)
    #endif
    } else{ // with B-TIS3 interpolation
        double *const __restrict__ BypartBTIS3  = smpi->dynamics_Bpart_yBTIS3[ithread].data();
        double *const __restrict__ BzpartBTIS3  = smpi->dynamics_Bpart_zBTIS3[ithread].data();
        const double *const __restrict__ By2D_mBTIS3 = static_cast<Field2D *>( EMfields->By_mBTIS3 )->data();
        const double *const __restrict__ Bz2D_mBTIS3 = static_cast<Field2D *>( EMfields->Bz_mBTIS3 )->data();
#if defined( SMILEI_ACCELERATOR_GPU_OMP )

    #pragma omp target map( to                                                     \
                            : i_domain_begin, j_domain_begin )                     \
        is_device_ptr /* map */ ( /* to: */                                        \
                                  position_x /* [first_index:npart_range_size] */, \
                                  position_y /* [first_index:npart_range_size] */ )
    #pragma omp teams distribute parallel for
#elif defined(SMILEI_ACCELERATOR_GPU_OACC)
    #pragma acc enter data create(this)
    #pragma acc update device(this)
    size_t interpolation_range_size = ( last_index + 1 * nparts ) - first_index;
    #pragma acc parallel present(ELoc [first_index:interpolation_range_size],\
                                 BLoc [first_index:interpolation_range_size],\
                                 BypartBTIS3 [first_index:interpolation_range_size],\
                                 BzpartBTIS3 [first_index:interpolation_range_size],\
                                 iold [first_index:interpolation_range_size],\
                                 delta [first_index:interpolation_range_size],\
                                 Ex2D [0:sizeofEx],\
                                 Ey2D [0:sizeofEy],\
                                 Ez2D [0:sizeofEz],\
                                 Bx2D [0:sizeofBx],\
                                 By2D [0:sizeofBy],\
                                 Bz2D [0:sizeofBz],\
                                 By2D_mBTIS3 [0:sizeofEz],\
                                 Bz2D_mBTIS3 [0:sizeofEy])\
    deviceptr(position_x, position_y)              \
    copyin(d_inv_[0:2])
    #pragma acc loop gang worker vector
#endif
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            // Normalized particle position
            const double xpn = position_x[ipart]*d_inv_[0];
            const double ypn = position_y[ipart]*d_inv_[1];

            // Calculate coeffs

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[3];
            double coeffxd[3], coeffyd[3];

            coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // Interpolation of Ex^(d,p)
            ELoc[0*nparts+ipart]          = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1], ny_p );
            // Interpolation of Ey^(p,d)
            ELoc[1*nparts+ipart]          = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1], ny_d );
            // Interpolation of Ez^(p,p)
            ELoc[2*nparts+ipart]          = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1], ny_p );
            // Interpolation of Bx^(p,d)
            BLoc[0*nparts+ipart]          = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1], ny_d );
            // Interpolation of By^(d,p)
            BLoc[1*nparts+ipart]          = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1], ny_p );
            // Interpolation of Bz^(d,d)
            BLoc[2*nparts+ipart]          = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1], ny_d );
            // Interpolation of ByBTIS3^(p,p)
            BypartBTIS3[0*nparts+ipart ]  = compute( &coeffxp[1], &coeffyp[1], By2D_mBTIS3, idx_p[0], idx_p[1], ny_p );
            // Interpolation of BzBTIS3^(p,d)
            BzpartBTIS3[0*nparts+ipart ]  = compute( &coeffxp[1], &coeffyd[1], Bz2D_mBTIS3, idx_p[0], idx_d[1], ny_d );

            //Buffering of iol and delta
            iold[0*nparts+ipart]  = idx_p[0];
            iold[1*nparts+ipart]  = idx_p[1];
            delta[0*nparts+ipart] = delta_p[0];
            delta[1*nparts+ipart] = delta_p[1];

        } // end ipart loop
    #if defined(SMILEI_ACCELERATOR_GPU_OACC)
        #pragma acc exit data delete(this)
    #endif
    } // end with B-TIS interpolation
}

// -----------------------------------------------------------------------------
//! Interpolator specific to tracked particles.
//! A selection of particles may be provided
// -----------------------------------------------------------------------------
void Interpolator2D2Order::fieldsSelection( ElectroMagn *EMfields,
                                            Particles &particles,
                                            double *buffer,
                                            int offset,
                                            vector<unsigned int> *selection )
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

void Interpolator2D2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the electromagnetic fields
    Field2D *Ex2D = static_cast<Field2D *>( EMfields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( EMfields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( EMfields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( EMfields->Bx_m );
    Field2D *By2D = static_cast<Field2D *>( EMfields->By_m );
    Field2D *Bz2D = static_cast<Field2D *>( EMfields->Bz_m );

    // Static cast of the envelope fields
    Field2D *Phi2D = static_cast<Field2D *>( EMfields->envelope->Phi_ );
    Field2D *GradPhix2D = static_cast<Field2D *>( EMfields->envelope->GradPhix_ );
    Field2D *GradPhiy2D = static_cast<Field2D *>( EMfields->envelope->GradPhiy_ );
    Field2D *GradPhiz2D = static_cast<Field2D *>( EMfields->envelope->GradPhiz_ );

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    
    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            //fieldsForTasks( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );

            // Normalized particle position
            double xpn = particles.position( 0, ipart )*d_inv_[0];
            double ypn = particles.position( 1, ipart )*d_inv_[1];

            // Calculate coeffs

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[3];
            double coeffxd[3], coeffyd[3];

            coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // Interpolation of Ex^(d,p)
            ( *Epart )[ipart+0*nparts] = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1] );

            // Interpolation of Ey^(p,d)
            ( *Epart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1] );

            // Interpolation of Ez^(p,p)
            ( *Epart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1] );

            // Interpolation of Bx^(p,d)
            ( *Bpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1] );

            // Interpolation of By^(d,p)
            ( *Bpart )[ipart+1*nparts] = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1] );

            // Interpolation of Bz^(d,d)
            ( *Bpart )[ipart+2*nparts] = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1] );

            // -------------------------
            // Interpolation of Phi^(p,p)
            // -------------------------
            ( *PHIpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], Phi2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhix^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhix2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhiy^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiy2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhiz^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiz2D, idx_p[0], idx_p[1] );

            //Buffering of iol and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];

        }  // end ipart loop
      
    } else { // with B-TIS3 interpolation
      
        Field2D *By2DBTIS3 = static_cast<Field2D *>( EMfields->By_mBTIS3 );
        Field2D *Bz2DBTIS3 = static_cast<Field2D *>( EMfields->Bz_mBTIS3 );
        std::vector<double> *BpartyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread] );
        std::vector<double> *BpartzBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread] );
        
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            //fieldsForTasks( EMfields, particles, ipart, nparts, &( *Epart )[ipart], &( *Bpart )[ipart] );

            // Normalized particle position
            double xpn = particles.position( 0, ipart )*d_inv_[0];
            double ypn = particles.position( 1, ipart )*d_inv_[1];

            // Calculate coeffs

            int idx_p[2], idx_d[2];
            double delta_p[2];
            double coeffxp[3], coeffyp[3];
            double coeffxd[3], coeffyd[3];

            coeffs( xpn, ypn, idx_p, idx_d, coeffxp, coeffyp, coeffxd, coeffyd, delta_p );

            // Interpolation of Ex^(d,p)
            ( *Epart )[ipart+0*nparts] = compute( &coeffxd[1], &coeffyp[1], Ex2D, idx_d[0], idx_p[1] );

            // Interpolation of Ey^(p,d)
            ( *Epart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyd[1], Ey2D, idx_p[0], idx_d[1] );

            // Interpolation of Ez^(p,p)
            ( *Epart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], Ez2D, idx_p[0], idx_p[1] );

            // Interpolation of Bx^(p,d)
            ( *Bpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], Bx2D, idx_p[0], idx_d[1] );

            // Interpolation of By^(d,p)
            ( *Bpart )[ipart+1*nparts] = compute( &coeffxd[1], &coeffyp[1], By2D, idx_d[0], idx_p[1] );

            // Interpolation of Bz^(d,d)
            ( *Bpart )[ipart+2*nparts] = compute( &coeffxd[1], &coeffyd[1], Bz2D, idx_d[0], idx_d[1] );
            
            // Interpolation of ByBTIS3^(p,p)
            ( *BpartyBTIS3 )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], By2DBTIS3, idx_p[0], idx_p[1] );
            
            // Interpolation of BzBTIS3^(p,d)
            ( *BpartzBTIS3 )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], Bz2DBTIS3, idx_p[0], idx_d[1] );


            // -------------------------
            // Interpolation of Phi^(p,p)
            // -------------------------
            ( *PHIpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], Phi2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhix^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhix2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhiy^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiy2D, idx_p[0], idx_p[1] );

            // -------------------------
            // Interpolation of GradPhiz^(p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiz2D, idx_p[0], idx_p[1] );

            //Buffering of iol and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];

      }  // end ipart loop
      
    } // end with B-TIS3 interpolation
    


} // END Interpolator2D2OrderForTasks


void Interpolator2D2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field2D *Phi_m2D = static_cast<Field2D *>( EMfields->envelope->Phi_m );

    Field2D *GradPhix_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhix_m );
    Field2D *GradPhiy_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhiy_m );
    Field2D *GradPhiz_m2D = static_cast<Field2D *>( EMfields->envelope->GradPhiz_m );

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

        // Calculate coeffs
        int idx_p[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        coeffs( xpn, ypn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

        // -------------------------
        // Interpolation of Phiold^(p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], Phi_m2D, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhixold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhix_m2D, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhiyold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiy_m2D, idx_p[0], idx_p[1] );

        // -------------------------
        // Interpolation of GradPhizold^(p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], GradPhiz_m2D, idx_p[0], idx_p[1] );

        //Buffering of iol and delta
        ( *iold )[ipart+0*nparts]  = idx_p[0];
        ( *iold )[ipart+1*nparts]  = idx_p[1];
        ( *delta )[ipart+0*nparts] = delta_p[0];
        ( *delta )[ipart+1*nparts] = delta_p[1];


    }

} // END Interpolator2D2OrderForTasks


void Interpolator2D2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
{
    // Static cast of the electromagnetic fields
    Field2D *Env_A_abs_2D = static_cast<Field2D *>( EMfields->Env_A_abs_ );
    Field2D *Env_Chi_2D = static_cast<Field2D *>( EMfields->Env_Chi_ );
    Field2D *Env_E_abs_2D = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    Field2D *Env_Ex_abs_2D = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];

    int idx_p[2];
    double delta_p[2];
    double coeffxp[3], coeffyp[3];
    
    coeffs( xpn, ypn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

    // -------------------------
    // Interpolation of Env_A_abs_^(p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], Env_A_abs_2D, idx_p[0], idx_p[1] );

    // -------------------------
    // Interpolation of Env_Chi_^(p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp[1], &coeffyp[1], Env_Chi_2D, idx_p[0], idx_p[1] );

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], Env_E_abs_2D, idx_p[0], idx_p[1] );

    // -------------------------
    // Interpolation of Env_Ex_abs_^(p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], Env_Ex_abs_2D, idx_p[0], idx_p[1] );

} // END Interpolator2D2Order

void Interpolator2D2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Static cast of the envelope fields
    Field2D *EnvEabs = static_cast<Field2D *>( EMfields->Env_E_abs_ );
    Field2D *EnvExabs = static_cast<Field2D *>( EMfields->Env_Ex_abs_ );

    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        
        // Calculate coeffs
        int idx_p[2];
        double delta_p[2];
        double coeffxp[3], coeffyp[3];
        coeffs( xpn, ypn, idx_p, NULL, coeffxp, coeffyp, NULL, NULL, delta_p );

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], EnvEabs, idx_p[0], idx_p[1] );
        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], EnvExabs, idx_p[0], idx_p[1] );

    }


} // END Interpolator2D2Order
