#include <iostream>
#include <cmath>

#include "Species.h"
#include "PusherBorisNR.h"
#include "Particles.h"

using namespace std;

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
    
    
    int nparts;
    if (vecto) {
        nparts = Epart->size()/3;
    } else {
        nparts = particles.size();
    }
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    double* position_x = particles.getPtrPosition(0);
    double* position_y = NULL;
    double* position_z = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
        }
    }
    
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    
    short *charge = particles.getPtrCharge();
    
    double charge_over_mass ;
    double umx, umy, umz;
    double upx, upy, upz;
    double alpha;
    double Tx, Ty, Tz;
    double T2;
    double Sx, Sy, Sz;
    
    for( int ipart=istart ; ipart<iend; ipart++ ) {
    
        charge_over_mass = (double)( charge[ipart] ) *one_over_mass_;
        alpha = charge_over_mass*dts2;
        
        // uminus = v + q/m * dt/2 * E
        umx = momentum_x[ipart] * one_over_mass_ + alpha * ( *( Ex+ipart-ipart_buffer_offset ) );
        umy = momentum_y[ipart] * one_over_mass_ + alpha * ( *( Ey+ipart-ipart_buffer_offset ) );
        umz = momentum_z[ipart] * one_over_mass_ + alpha * ( *( Ez+ipart-ipart_buffer_offset ) );
        
        
        // Rotation in the magnetic field
        
        Tx    = alpha * ( *( Bx+ipart-ipart_buffer_offset ) );
        Ty    = alpha * ( *( By+ipart-ipart_buffer_offset ) );
        Tz    = alpha * ( *( Bz+ipart-ipart_buffer_offset ) );
        
        T2 = Tx*Tx + Ty*Ty + Tz*Tz;
        
        Sx = 2*Tx/( 1.+T2 );
        Sy = 2*Ty/( 1.+T2 );
        Sz = 2*Tz/( 1.+T2 );
        
        // uplus = uminus + uprims x S
        upx = umx + umy*Sz - umz*Sy;
        upy = umy + umz*Sx - umx*Sz;
        upz = umz + umx*Sy - umy*Sx;
        
        
        momentum_x[ipart] = mass_ * ( upx + alpha*( *( Ex+ipart-ipart_buffer_offset ) ) );
        momentum_y[ipart] = mass_ * ( upy + alpha*( *( Ey+ipart-ipart_buffer_offset ) ) );
        momentum_z[ipart] = mass_ * ( upz + alpha*( *( Ez+ipart-ipart_buffer_offset ) ) );
        
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
