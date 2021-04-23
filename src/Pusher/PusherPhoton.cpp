// -----------------------------------------------------------------------------
//! \file PusherHigueraCary.cpp
//
//!  \brief Class fot the photon pusher
//
//  @date 2017-07-24
// -----------------------------------------------------------------------------

#include "PusherPhoton.h"

#include <iostream>
#include <cmath>

#include "Species.h"
#include "Particles.h"

using namespace std;

PusherPhoton::PusherPhoton( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherPhoton::~PusherPhoton()
{
}

/***********************************************************************
    Rectilinear propagation of the photons
***********************************************************************/

void PusherPhoton::operator()( Particles &particles, SmileiMPI *smpi,
                               int istart, int iend, int ithread, int ipart_ref )
{
    // Inverse normalized energy
    std::vector<double> *invgf = &( smpi->dynamics_invgf[ithread] );
    
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, 0 ) );
    }
    double *position[3];
    for( int i = 0 ; i<nDim_ ; i++ ) {
        position[i] =  &( particles.position( i, 0 ) );
    }
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
    
        ( *invgf )[ipart] = 1. / sqrt( momentum[0][ipart]*momentum[0][ipart] +
                                       momentum[1][ipart]*momentum[1][ipart] +
                                       momentum[2][ipart]*momentum[2][ipart] );
                                       
        // Move the photons
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i][ipart]     += dt*momentum[i][ipart]*( *invgf )[ipart];
        }
        
    }
    
    if( vecto ) {
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
