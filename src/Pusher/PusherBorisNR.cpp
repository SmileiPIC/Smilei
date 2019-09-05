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

void PusherBorisNR::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    
    int nparts = particles.size();
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    double charge_over_mass_ ;
    double umx, umy, umz;
    double upx, upy, upz;
    double alpha;
    double Tx, Ty, Tz;
    double T2;
    double Sx, Sy, Sz;
    
    for( int ipart=istart ; ipart<iend; ipart++ ) {
    
        charge_over_mass_ = static_cast<double>( particles.charge( ipart ) )*one_over_mass_;
        alpha = charge_over_mass_*dts2;
        
        // uminus = v + q/m * dt/2 * E
        umx = particles.momentum( 0, ipart ) * one_over_mass_ + alpha * ( *( Ex+ipart ) );
        umy = particles.momentum( 1, ipart ) * one_over_mass_ + alpha * ( *( Ey+ipart ) );
        umz = particles.momentum( 2, ipart ) * one_over_mass_ + alpha * ( *( Ez+ipart ) );
        
        
        // Rotation in the magnetic field
        
        Tx    = alpha * ( *( Bx+ipart ) );
        Ty    = alpha * ( *( By+ipart ) );
        Tz    = alpha * ( *( Bz+ipart ) );
        
        T2 = Tx*Tx + Ty*Ty + Tz*Tz;
        
        Sx = 2*Tx/( 1.+T2 );
        Sy = 2*Ty/( 1.+T2 );
        Sz = 2*Tz/( 1.+T2 );
        
        // uplus = uminus + uprims x S
        upx = umx + umy*Sz - umz*Sy;
        upy = umy + umz*Sx - umx*Sz;
        upz = umz + umx*Sy - umy*Sx;
        
        
        particles.momentum( 0, ipart ) = mass_ * ( upx + alpha*( *( Ex+ipart ) ) );
        particles.momentum( 1, ipart ) = mass_ * ( upy + alpha*( *( Ey+ipart ) ) );
        particles.momentum( 2, ipart ) = mass_ * ( upz + alpha*( *( Ez+ipart ) ) );
        
        // Move the particle
        for( int i = 0 ; i<nDim_ ; i++ ) {
            particles.position( i, ipart )     += dt*particles.momentum( i, ipart );
        }
    }
    
    if( vecto ) {
        double *position[3];
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i] =  &( particles.position( i, 0 ) );
        }
        int *cell_keys;
        particles.cell_keys.resize( iend-istart );
        cell_keys = &( particles.cell_keys[0] );
        
        #pragma omp simd
        for( int ipart=istart ; ipart<iend; ipart++ ) {
        
            for( int i = 0 ; i<nDim_ ; i++ ) {
                cell_keys[ipart] *= nspace[i];
                cell_keys[ipart] += round( ( position[i][ipart]-min_loc_vec[i] ) * dx_inv_[i] );
            }
            
        }
    }
    
}
