#include "Interpolator3D2Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field3D.h"
#include "LaserEnvelope.h"
#include "Particles.h"


// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator3D2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator3D2Order::Interpolator3D2Order( Params &params, Patch *patch ) : Interpolator3D( patch )
{

    d_inv_[0] = 1.0/params.cell_length[0];
    d_inv_[1] = 1.0/params.cell_length[1];
    d_inv_[2] = 1.0/params.cell_length[2];

}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator3D2Order::fields( ElectroMagn *EMfields, Particles &particles, int ipart, int nparts, double *ELoc, double *BLoc )
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
    double coeffxp[3], coeffyp[3], coeffzp[3];
    double coeffxd[3], coeffyd[3], coeffzd[3];
    coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Ey3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2] );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2] );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2] );
} // END Interpolator3D2Order

void Interpolator3D2Order::fieldsAndCurrents( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *, int ithread, LocalFields *JLoc, double *RhoLoc )
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
    Field3D *By3DBTIS3;
    Field3D *Bz3DBTIS3;
    if(smpi->use_BTIS3){
        By3DBTIS3 = static_cast<Field3D *>( EMfields->By_mBTIS3 );
        Bz3DBTIS3 = static_cast<Field3D *>( EMfields->Bz_mBTIS3 );
    }

    // Normalized particle position
    double xpn = particles.position( 0, ipart )*d_inv_[0];
    double ypn = particles.position( 1, ipart )*d_inv_[1];
    double zpn = particles.position( 2, ipart )*d_inv_[2];
    // Compute coeffs
    int    idx_p[3], idx_d[3];
    double delta_p[3];
    double coeffxp[3], coeffyp[3], coeffzp[3];
    double coeffxd[3], coeffyd[3], coeffzd[3];
    coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

    int nparts( particles.numberOfParticles() );

    // Interpolation of Ex^(d,p,p)
    *( ELoc+0*nparts ) = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Ey^(p,d,p)
    *( ELoc+1*nparts ) = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Ey3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Ez^(p,p,d)
    *( ELoc+2*nparts ) = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Bx^(p,d,d)
    *( BLoc+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2] );
    // Interpolation of By^(d,p,d)
    *( BLoc+1*nparts ) = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2] );
    // Interpolation of Bz^(d,d,p)
    *( BLoc+2*nparts ) = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2] );
    // Interpolation of Jx^(d,p,p)
    JLoc->x = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Jx3D, idx_d[0], idx_p[1], idx_p[2] );
    // Interpolation of Jy^(p,d,p)
    JLoc->y = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Jy3D, idx_p[0], idx_d[1], idx_p[2] );
    // Interpolation of Jz^(p,p,d)
    JLoc->z = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Jz3D, idx_p[0], idx_p[1], idx_d[2] );
    // Interpolation of Rho^(p,p,p)
    ( *RhoLoc ) = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Rho3D, idx_p[0], idx_p[1], idx_p[2] );
    
    if (smpi->use_BTIS3){
        // Interpolation of ByBTIS3^(p,p,d)
        *( BLocyBTIS3+0*nparts ) = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], By3DBTIS3, idx_p[0], idx_p[1], idx_d[2] );
        // Interpolation of BzBTIS3^(p,d,p)
        *( BLoczBTIS3+0*nparts ) = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Bz3DBTIS3, idx_p[0], idx_d[1], idx_p[2] );
    }

}

// Interpolator on another field than the basic ones
void Interpolator3D2Order::oneField( Field **field, Particles &particles, int *istart, int *iend, double *FieldLoc, double *, double *, double * )
{
    Field3D *F = static_cast<Field3D *>( *field );
    int    idx_p[3], idx_d[3];
    double delta_p[3];
    double coeffxp[3], coeffyp[3], coeffzp[3];
    double coeffxd[3], coeffyd[3], coeffzd[3];
    double *coeffx = F->isDual( 0 ) ? &coeffxd[1] : &coeffxp[1];
    double *coeffy = F->isDual( 1 ) ? &coeffyd[1] : &coeffyp[1];
    double *coeffz = F->isDual( 2 ) ? &coeffzd[1] : &coeffzp[1];
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

void Interpolator3D2Order::fieldsWrapper( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, unsigned int, int )
{
    const int nparts = particles.numberOfParticles();

    // CCE 13 implementation of OpenMP (as of 2022/04/07) does not like
    // dereferenced ptrs in the for loop's condition.
    const int first_index = *istart;
    const int last_index  = *iend;

    if( first_index == last_index ) {
        // Early exit to avoid UB. If the dynamics_* array are empty() we can't
        // call data().
        return;
    }

    SMILEI_ASSERT( !smpi->dynamics_Epart[ithread].empty() );
    SMILEI_ASSERT( !smpi->dynamics_Bpart[ithread].empty() );
    SMILEI_ASSERT( !smpi->dynamics_iold[ithread].empty() );
    SMILEI_ASSERT( !smpi->dynamics_deltaold[ithread].empty() );

    double *const __restrict__ ELoc  = smpi->dynamics_Epart[ithread].data();
    double *const __restrict__ BLoc  = smpi->dynamics_Bpart[ithread].data();
    int *const __restrict__ iold     = smpi->dynamics_iold[ithread].data();
    double *const __restrict__ delta = smpi->dynamics_deltaold[ithread].data();

    const double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    const double *const __restrict__ position_y = particles.getPtrPosition( 1 );
    const double *const __restrict__ position_z = particles.getPtrPosition( 2 );

    const double *const __restrict__ Ex3D = EMfields->Ex_->data_;
    const double *const __restrict__ Ey3D = EMfields->Ey_->data_;
    const double *const __restrict__ Ez3D = EMfields->Ez_->data_;
    const double *const __restrict__ Bx3D = EMfields->Bx_m->data_;
    const double *const __restrict__ By3D = EMfields->By_m->data_;
    const double *const __restrict__ Bz3D = EMfields->Bz_m->data_;

#if defined(SMILEI_ACCELERATOR_GPU_OACC)
    const int sizeofEx = EMfields->Ex_->size();
    const int sizeofEy = EMfields->Ey_->size();
    const int sizeofEz = EMfields->Ez_->size();
    const int sizeofBx = EMfields->Bx_m->size();
    const int sizeofBy = EMfields->By_m->size();
    const int sizeofBz = EMfields->Bz_m->size();
#endif

    const int nx_p = EMfields->Bx_m->dims_[0];
    const int ny_p = EMfields->By_m->dims_[1];
    const int nz_p = EMfields->Bz_m->dims_[2];
    const int nx_d = nx_p + 1;
    const int ny_d = ny_p + 1;
    const int nz_d = nz_p + 1;

    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
#if defined(SMILEI_ACCELERATOR_GPU_OMP)

    #pragma omp target map( to                                                 \
                            : i_domain_begin, j_domain_begin, k_domain_begin ) \
        is_device_ptr( /* to: */                                               \
                       position_x /* [first_index:npart_range_size] */,        \
                       position_y /* [first_index:npart_range_size] */,        \
                       position_z /* [first_index:npart_range_size] */ )
    #pragma omp teams distribute parallel for
#elif defined(SMILEI_ACCELERATOR_GPU_OACC)
    #pragma acc enter data create(this)
    #pragma acc update device(this)
    size_t interpolation_range_size = ( last_index + 2 * nparts ) - first_index;
    #pragma acc parallel present(ELoc [first_index:interpolation_range_size],  \
                                 BLoc [first_index:interpolation_range_size],  \
                                 iold [first_index:interpolation_range_size],  \
                                delta [first_index:interpolation_range_size], \
                                 Ex3D [0:sizeofEx],                            \
                                 Ey3D [0:sizeofEy],                            \
                                 Ez3D [0:sizeofEz],                            \
                                 Bx3D [0:sizeofBx],                            \
                                 By3D [0:sizeofBy],                            \
                                 Bz3D [0:sizeofBz])                            \
        deviceptr(position_x,                                                  \
                  position_y,                                                  \
                  position_z)                                                  \
        copyin(d_inv_[0:3])

    #pragma acc loop gang worker vector
#endif
        for( int ipart=first_index ; ipart<last_index; ipart++ ) {

            // Interpolation on current particle

            // Normalized particle position
            const double xpn = position_x[ipart]*d_inv_[0];
            const double ypn = position_y[ipart]*d_inv_[1];
            const double zpn = position_z[ipart]*d_inv_[2];

            // Compute coeffs
            int    idx_p[3], idx_d[3];
            double delta_p[3];
            double coeffxp[3], coeffyp[3], coeffzp[3];
            double coeffxd[3], coeffyd[3], coeffzd[3];
            coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

            // Interpolation of Ex^(d,p,p)
            ELoc[0*nparts+ipart] = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2], nx_d, ny_p, nz_p );
            // Interpolation of Ey^(p,d,p)
            ELoc[1*nparts+ipart] = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Ey3D, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );
            // Interpolation of Ez^(p,p,d)
            ELoc[2*nparts+ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );
            // Interpolation of Bx^(p,d,d)
            BLoc[0*nparts+ipart] = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2], nx_p, ny_d, nz_d );
            // Interpolation of By^(d,p,d)
            BLoc[1*nparts+ipart] = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2], nx_d, ny_p, nz_d );
            // Interpolation of Bz^(d,d,p)
            BLoc[2*nparts+ipart] = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2], nx_d, ny_d, nz_p );

            // Buffering of iol and delta
            iold[0*nparts+ipart]  = idx_p[0];
            iold[1*nparts+ipart]  = idx_p[1];
            iold[2*nparts+ipart]  = idx_p[2];
            delta[0*nparts+ipart] = delta_p[0];
            delta[1*nparts+ipart] = delta_p[1];
            delta[2*nparts+ipart] = delta_p[2];
        }
        #if defined(SMILEI_ACCELERATOR_GPU_OACC)
            #pragma acc exit data delete(this)
        #endif
    } else { // with B-TIS3 interpolation

        const double *const __restrict__ By3DBTIS3 = EMfields->By_mBTIS3->data_;
        const double *const __restrict__ Bz3DBTIS3 = EMfields->Bz_mBTIS3->data_;
        double *const __restrict__ BLocyBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread][0] );
        double *const __restrict__ BLoczBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread][0] );
        
        // loop on particles from istart to iend
#if defined(SMILEI_ACCELERATOR_GPU_OMP)

    #pragma omp target map( to                                                 \
                            : i_domain_begin, j_domain_begin, k_domain_begin ) \
        is_device_ptr( /* to: */                                               \
                       position_x /* [first_index:npart_range_size] */,        \
                       position_y /* [first_index:npart_range_size] */,        \
                       position_z /* [first_index:npart_range_size] */ )
    #pragma omp teams distribute parallel for
#elif defined(SMILEI_ACCELERATOR_GPU_OACC)
    #pragma acc enter data create(this)
    #pragma acc update device(this)
    size_t interpolation_range_size = ( last_index + 2 * nparts ) - first_index;
    #pragma acc parallel present(ELoc [first_index:interpolation_range_size],  \
                                 BLoc [first_index:interpolation_range_size],  \
                                 BLocyBTIS3 [first_index:interpolation_range_size],\
                                 BLoczBTIS3 [first_index:interpolation_range_size],\
                                 iold [first_index:interpolation_range_size],  \
                                delta [first_index:interpolation_range_size], \
                                 Ex3D [0:sizeofEx],                            \
                                 Ey3D [0:sizeofEy],                            \
                                 Ez3D [0:sizeofEz],                            \
                                 Bx3D [0:sizeofBx],                            \
                                 By3D [0:sizeofBy],                            \
                                 Bz3D [0:sizeofBz],                            \
                                 By3DBTIS3 [0:sizeofEz],                     \
                                 Bz3DBTIS3 [0:sizeofEy])                     \
        deviceptr(position_x,                                                  \
                  position_y,                                                  \
                  position_z)                                                  \
        copyin(d_inv_[0:3])

    #pragma acc loop gang worker vector
#endif
        for( int ipart=first_index ; ipart<last_index; ipart++ ) {

            //Interpolation on current particle

            // Normalized particle position
            const double xpn = position_x[ ipart ]*d_inv_[0];
            const double ypn = position_y[ ipart ]*d_inv_[1];
            const double zpn = position_z[ ipart ]*d_inv_[2];
            // Compute coeffs

            int idx_p[3], idx_d[3];
            double delta_p[3];
            double coeffxp[3], coeffyp[3], coeffzp[3];
            double coeffxd[3], coeffyd[3], coeffzd[3];

            coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

            // Interpolation of Ex^(d,p,p)
            *( ELoc+0*nparts+ipart ) = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2], nx_d, ny_p, nz_p );
            // Interpolation of Ey^(p,d,p)
            *( ELoc+1*nparts+ipart ) = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Ey3D, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );
            // Interpolation of Ez^(p,p,d)
            *( ELoc+2*nparts+ipart ) = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );
            // Interpolation of Bx^(p,d,d)
            *( BLoc+0*nparts+ipart ) = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2], nx_p, ny_d, nz_d );
            // Interpolation of By^(d,p,d)
            *( BLoc+1*nparts+ipart ) = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2], nx_d, ny_p, nz_d );
            // Interpolation of Bz^(d,d,p)
            *( BLoc+2*nparts+ipart ) = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2], nx_d, ny_d, nz_p );
            // Interpolation of By^(p,p,d)
            *( BLocyBTIS3+0*nparts+ipart ) = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], By3DBTIS3, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );
            // Interpolation of By^(p,d,p)
            *( BLoczBTIS3+0*nparts+ipart ) = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Bz3DBTIS3, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );

            //Buffering of iol and delta
            iold[ipart+0*nparts]  = idx_p[0];
            iold[ipart+1*nparts]  = idx_p[1];
            iold[ipart+2*nparts]  = idx_p[2];
            delta[ipart+0*nparts] = delta_p[0];
            delta[ipart+1*nparts] = delta_p[1];
            delta[ipart+2*nparts] = delta_p[2];
            #if defined(SMILEI_ACCELERATOR_GPU_OACC)
                #pragma acc exit data delete(this)
            #endif
        } // end ipart loop
      
    } // end with B-TIS3 interpolation
}


// Interpolator specific to tracked particles. A selection of particles may be provided
void Interpolator3D2Order::fieldsSelection( ElectroMagn *EMfields, Particles &particles, double *buffer, int offset, std::vector<unsigned int> *selection )
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

void Interpolator3D2Order::fieldsAndEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Electromagnetic fields
    double* Ex3D = EMfields->Ex_->data_;
    double* Ey3D = EMfields->Ey_->data_;
    double* Ez3D = EMfields->Ez_->data_;
    double* Bx3D = EMfields->Bx_m->data_;
    double* By3D = EMfields->By_m->data_;
    double* Bz3D = EMfields->Bz_m->data_;

    // Envelope fields
    double* Phi3D = EMfields->envelope->Phi_->data_;
    double* GradPhix3D = EMfields->envelope->GradPhix_->data_;
    double* GradPhiy3D = EMfields->envelope->GradPhiy_->data_;
    double* GradPhiz3D = EMfields->envelope->GradPhiz_->data_;

    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    std::vector<double> *PHIpart        = &( smpi->dynamics_PHIpart[ithread] );
    std::vector<double> *GradPHIpart    = &( smpi->dynamics_GradPHIpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    int nx_p = EMfields->Bx_m->dims_[0];
    int ny_p = EMfields->By_m->dims_[1];
    int nz_p = EMfields->Bz_m->dims_[2];
    int nx_d = nx_p+1;
    int ny_d = ny_p+1;
    int nz_d = nz_p+1;

    //Loop on bin particles
    int nparts( particles.numberOfParticles() );
    if (!smpi->use_BTIS3){ // without B-TIS3 interpolation
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            int idx_p[3], idx_d[3];
            double delta_p[3];
            double coeffxp[3], coeffyp[3], coeffzp[3];
            double coeffxd[3], coeffyd[3], coeffzd[3];

            // Normalized particle position
            double xpn = particles.position( 0, ipart )*d_inv_[0];
            double ypn = particles.position( 1, ipart )*d_inv_[1];
            double zpn = particles.position( 2, ipart )*d_inv_[2];
            coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

            // Interpolation of Ex^(d,p,p)
            ( *Epart )[ipart+0*nparts] = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2], nx_d, ny_p, nz_p );

            // Interpolation of Ey^(p,d,p)
            ( *Epart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1],  Ey3D, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );

            // Interpolation of Ez^(p,p,d)
            ( *Epart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );

            // Interpolation of Bx^(p,d,d)
            ( *Bpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2], nx_p, ny_d, nz_d );

            // Interpolation of By^(d,p,d)
            ( *Bpart )[ipart+1*nparts] = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2], nx_d, ny_p, nz_d );

            // Interpolation of Bz^(d,d,p)
            ( *Bpart )[ipart+2*nparts] = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2], nx_d, ny_d, nz_p );


            // -------------------------
            // Interpolation of Phi^(p,p,p)
            // -------------------------
            ( *PHIpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Phi3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhix^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhix3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhiy^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiy3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhiz^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiz3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            //Buffering of iol and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *iold )[ipart+2*nparts]  = idx_p[2];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];
            ( *delta )[ipart+2*nparts] = delta_p[2];

        }  // end ipart loop
      
    } else { // with B-TIS3 interpolation
      
        std::vector<double> *BypartBTIS3;
        std::vector<double> *BzpartBTIS3;

        BypartBTIS3 = &( smpi->dynamics_Bpart_yBTIS3[ithread] );
        BzpartBTIS3 = &( smpi->dynamics_Bpart_zBTIS3[ithread] );
        double* By3DBTIS3 = EMfields->By_mBTIS3->data_;
        double* Bz3DBTIS3 = EMfields->Bz_mBTIS3->data_;
      
        for( int ipart=*istart ; ipart<*iend; ipart++ ) {

            int idx_p[3], idx_d[3];
            double delta_p[3];
            double coeffxp[3], coeffyp[3], coeffzp[3];
            double coeffxd[3], coeffyd[3], coeffzd[3];

            // Normalized particle position
            double xpn = particles.position( 0, ipart )*d_inv_[0];
            double ypn = particles.position( 1, ipart )*d_inv_[1];
            double zpn = particles.position( 2, ipart )*d_inv_[2];

            coeffs( xpn, ypn, zpn, idx_p, idx_d, coeffxp, coeffyp, coeffzp, coeffxd, coeffyd, coeffzd, delta_p );

            // Interpolation of Ex^(d,p,p)
            ( *Epart )[ipart+0*nparts] = compute( &coeffxd[1], &coeffyp[1], &coeffzp[1], Ex3D, idx_d[0], idx_p[1], idx_p[2], nx_d, ny_p, nz_p );

            // Interpolation of Ey^(p,d,p)
            ( *Epart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1],  Ey3D, idx_p[0], idx_d[1], idx_p[2], nx_p, ny_d, nz_p );

            // Interpolation of Ez^(p,p,d)
            ( *Epart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], Ez3D, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );

            // Interpolation of Bx^(p,d,d)
            ( *Bpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], &coeffzd[1], Bx3D, idx_p[0], idx_d[1], idx_d[2], nx_p, ny_d, nz_d );

            // Interpolation of By^(d,p,d)
            ( *Bpart )[ipart+1*nparts] = compute( &coeffxd[1], &coeffyp[1], &coeffzd[1], By3D, idx_d[0], idx_p[1], idx_d[2], nx_d, ny_p, nz_d );

            // Interpolation of Bz^(d,d,p)
            ( *Bpart )[ipart+2*nparts] = compute( &coeffxd[1], &coeffyd[1], &coeffzp[1], Bz3D, idx_d[0], idx_d[1], idx_p[2], nx_d, ny_d, nz_p );
            
            // Interpolation of ByBTIS3^(p,p,d)
            (*BypartBTIS3)[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzd[1], By3DBTIS3, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d );
            
            // Interpolation of BzBTIS3^(p,d,p)
            (*BzpartBTIS3)[ipart+0*nparts] = compute( &coeffxp[1], &coeffyd[1], &coeffzp[1], Bz3DBTIS3, idx_p[0], idx_p[1], idx_d[2], nx_p, ny_p, nz_d);


            // -------------------------
            // Interpolation of Phi^(p,p,p)
            // -------------------------
            ( *PHIpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Phi3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhix^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhix3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhiy^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiy3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );

            // -------------------------
            // Interpolation of GradPhiz^(p,p,p)
            // -------------------------
            ( *GradPHIpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiz3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p );


            //Buffering of iol and delta
            ( *iold )[ipart+0*nparts]  = idx_p[0];
            ( *iold )[ipart+1*nparts]  = idx_p[1];
            ( *iold )[ipart+2*nparts]  = idx_p[2];
            ( *delta )[ipart+0*nparts] = delta_p[0];
            ( *delta )[ipart+1*nparts] = delta_p[1];
            ( *delta )[ipart+2*nparts] = delta_p[2];

        }  // end ipart loop 
      
    } // end withB-TIS3 interpolation 
    
} // END Interpolator3D2Order

void Interpolator3D2Order::timeCenteredEnvelope( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Envelope fields
    double* Phi_m3D = EMfields->envelope->Phi_m->data_;
    double* GradPhix_m3D = EMfields->envelope->GradPhix_m->data_;
    double* GradPhiy_m3D = EMfields->envelope->GradPhiy_m->data_;
    double* GradPhiz_m3D = EMfields->envelope->GradPhiz_m->data_;

    std::vector<double> *PHI_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPHI_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );

    std::vector<int>    *iold  = &( smpi->dynamics_iold[ithread] );
    std::vector<double> *delta = &( smpi->dynamics_deltaold[ithread] );

    int nx_p = EMfields->Bx_m->dims_[0];
    int ny_p = EMfields->By_m->dims_[1];
    int nz_p = EMfields->Bz_m->dims_[2];

    //Loop on bin particles
    int nparts( particles.numberOfParticles());
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];

        int idx_p[3];
        double delta_p[3];
        double coeffxp[3], coeffyp[3], coeffzp[3];
        coeffs( xpn, ypn, zpn, idx_p, NULL, coeffxp, coeffyp, coeffzp, NULL, NULL, NULL, delta_p );

        // -------------------------
        // Interpolation of Phi_m^(p,p,p)
        // -------------------------
        ( *PHI_mpart )[ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Phi_m3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );

        // -------------------------
        // Interpolation of GradPhix_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+0*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhix_m3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );

        // -------------------------
        // Interpolation of GradPhiy_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+1*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiy_m3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );

        // -------------------------
        // Interpolation of GradPhiz_m^(p,p,p)
        // -------------------------
        ( *GradPHI_mpart )[ipart+2*nparts] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], GradPhiz_m3D, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );

        //Buffering of iold and delta
        ( *iold )[ipart+0*nparts]  = idx_p[0];
        ( *iold )[ipart+1*nparts]  = idx_p[1];
        ( *iold )[ipart+2*nparts]  = idx_p[2];
        ( *delta )[ipart+0*nparts] = delta_p[0];
        ( *delta )[ipart+1*nparts] = delta_p[1];
        ( *delta )[ipart+2*nparts] = delta_p[2];

    }

} // END Interpolator3D2Order

void Interpolator3D2Order::envelopeAndSusceptibility( ElectroMagn *EMfields, Particles &particles, int ipart, double *Env_A_abs_Loc, double *Env_Chi_Loc, double *Env_E_abs_Loc, double *Env_Ex_abs_Loc )
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

    int idx_p[3];
    double delta_p[3];
    double coeffxp[3], coeffyp[3], coeffzp[3];
    coeffs( xpn, ypn, zpn, idx_p, NULL, coeffxp, coeffyp, coeffzp, NULL, NULL, NULL, delta_p );

    // -------------------------
    // Interpolation of Env_A_abs_^(p,p,p)
    // -------------------------
    *( Env_A_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Env_A_abs_3D, idx_p[0], idx_p[1], idx_p[2] );

    // -------------------------
    // Interpolation of Env_Chi_^(p,p,p)
    // -------------------------
    *( Env_Chi_Loc ) = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Env_Chi_3D, idx_p[0], idx_p[1], idx_p[2] );

    // -------------------------
    // Interpolation of Env_E_abs_^(p,p,p)
    // -------------------------
    *( Env_E_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Env_E_abs_3D, idx_p[0], idx_p[1], idx_p[2] );


    // -------------------------
    // Interpolation of Env_E_abs_^(p,p,p)
    // -------------------------
    *( Env_Ex_abs_Loc ) = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], Env_Ex_abs_3D, idx_p[0], idx_p[1], idx_p[2] );

} // END Interpolator3D2Order

void Interpolator3D2Order::envelopeFieldForIonization( ElectroMagn *EMfields, Particles &particles, SmileiMPI *smpi, int *istart, int *iend, int ithread, int )
{
    // Envelope fields
    double* EnvEabs  = EMfields->Env_E_abs_->data_;
    double* EnvExabs = EMfields->Env_Ex_abs_->data_;

    std::vector<double> *EnvEabs_part  = &( smpi->dynamics_EnvEabs_part[ithread] );
    std::vector<double> *EnvExabs_part = &( smpi->dynamics_EnvExabs_part[ithread] );

    int nx_p = EMfields->Bx_m->dims_[0];
    int ny_p = EMfields->By_m->dims_[1];
    int nz_p = EMfields->Bz_m->dims_[2];

    //Loop on bin particles
    for( int ipart=*istart ; ipart<*iend; ipart++ ) {

        // Normalized particle position
        double xpn = particles.position( 0, ipart )*d_inv_[0];
        double ypn = particles.position( 1, ipart )*d_inv_[1];
        double zpn = particles.position( 2, ipart )*d_inv_[2];

        int idx_p[3];
        double delta_p[3];
        double coeffxp[3], coeffyp[3], coeffzp[3];
        coeffs( xpn, ypn, zpn, idx_p, NULL, coeffxp, coeffyp, coeffzp, NULL, NULL, NULL, delta_p );

        // ---------------------------------
        // Interpolation of Env_E_abs^(p,p,p)
        // ---------------------------------
        ( *EnvEabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], EnvEabs, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );

        // ---------------------------------
        // Interpolation of Env_Ex_abs^(p,p,p)
        // ---------------------------------
        ( *EnvExabs_part )[ipart] = compute( &coeffxp[1], &coeffyp[1], &coeffzp[1], EnvExabs, idx_p[0], idx_p[1], idx_p[2], nx_p, ny_p, nz_p  );
    }

} // END Interpolator3D2Order
