#include <iostream>
#include <cmath>

#include "Species.h"
#include "PusherBorisNR.h"
#include "Particles.h"

PusherBorisNR::PusherBorisNR( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherBorisNR::~PusherBorisNR()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBorisNR::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );

    const int nparts = vecto ? Epart->size() / 3 :
                               particles.last_index.back(); // particles.size()

    const double *const __restrict__ Ex = &( ( *Epart )[0*nparts] );
    const double *const __restrict__ Ey = &( ( *Epart )[1*nparts] );
    const double *const __restrict__ Ez = &( ( *Epart )[2*nparts] );
    const double *const __restrict__ Bx = &( ( *Bpart )[0*nparts] );
    const double *const __restrict__ By = &( ( *Bpart )[1*nparts] );
    const double *const __restrict__ Bz = &( ( *Bpart )[2*nparts] );

    double *const __restrict__ position_x = particles.getPtrPosition( 0 );
    double *const __restrict__ position_y = nDim_ > 1 ? particles.getPtrPosition( 1 ) : nullptr;
    double *const __restrict__ position_z = nDim_ > 2 ? particles.getPtrPosition( 2 ) : nullptr;

    double *const __restrict__ momentum_x = particles.getPtrMomentum( 0 );
    double *const __restrict__ momentum_y = particles.getPtrMomentum( 1 );
    double *const __restrict__ momentum_z = particles.getPtrMomentum( 2 );

    const short *__restrict__ charge = particles.getPtrCharge();

    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {

        const int ipart2 = ipart - ipart_buffer_offset;

        const double charge_over_mass = ( double )( charge[ipart] )*one_over_mass_;
        
        double alpha = charge_over_mass*dts2;

        // uminus = v + q/m * dt/2 * E
        const double umx = momentum_x[ipart] * one_over_mass_ + alpha * ( Ex[ipart2] );
        const double umy = momentum_y[ipart] * one_over_mass_ + alpha * ( Ey[ipart2] );
        const double umz = momentum_z[ipart] * one_over_mass_ + alpha * ( Ez[ipart2] );

        // Rotation in the magnetic field

        const double Tx    = alpha * ( Bx[ipart2] );
        const double Ty    = alpha * ( By[ipart2] );
        const double Tz    = alpha * ( Bz[ipart2] );

        const double T2 = Tx*Tx + Ty*Ty + Tz*Tz;

        const double Sx = 2*Tx/( 1.+T2 );
        const double Sy = 2*Ty/( 1.+T2 );
        const double Sz = 2*Tz/( 1.+T2 );

        // uplus = uminus + uprims x S
        const double upx = umx + umy*Sz - umz*Sy;
        const double upy = umy + umz*Sx - umx*Sz;
        const double upz = umz + umx*Sy - umy*Sx;


        momentum_x[ipart] = mass_ * ( upx + alpha*( Ex[ipart2] ) );
        momentum_y[ipart] = mass_ * ( upy + alpha*( Ey[ipart2] ) );
        momentum_z[ipart] = mass_ * ( upz + alpha*( Ez[ipart2] ) );

        // Move the particle
        position_x[ipart] += dt * momentum_x[ipart];
        if (nDim_>1) {
            position_y[ipart] += dt * momentum_y[ipart];
            if (nDim_>2) {
                position_z[ipart] += dt * momentum_z[ipart];
            }
        }

    }
    //
    // if( vecto ) {
    //     double *position[3];
    //     for( int i = 0 ; i<nDim_ ; i++ ) {
    //         position[i] =  &( particles.position( i, 0 ) );
    //     }
    //     int *cell_keys;
    //     particles.cell_keys.resize( iend-istart );
    //     cell_keys = &( particles.cell_keys[0] );
    //
    //     #pragma omp simd
    //     for( int ipart=istart ; ipart<iend; ipart++ ) {
    //
    //         for( int i = 0 ; i<nDim_ ; i++ ) {
    //             cell_keys[ipart] *= nspace[i];
    //             cell_keys[ipart] += round( ( position[i][ipart]-min_loc_vec[i] ) * dx_inv_[i] );
    //         }
    //
    //     }
    // }

}
