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

    PusherPhoton::PusherPhoton(Params& params, Species *species)
: Pusher(params, species)
{
}

PusherPhoton::~PusherPhoton()
{
}

/***********************************************************************
    Rectilinear propagation of the photons
***********************************************************************/

void PusherPhoton::operator() (Particles &particles, SmileiMPI* smpi,
                              int istart, int iend, int ithread)
{
    // Inverse normalized energy
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );
    double* position[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position[i] =  &( particles.position(i,0) );
#ifdef  __DEBUG
    double* position_old[3];
    for ( int i = 0 ; i<nDim_ ; i++ )
        position_old[i] =  &( particles.position_old(i,0) );
#endif

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) {

        (*invgf)[ipart] = 1. / sqrt( momentum[0][ipart]*momentum[0][ipart] +
                                     momentum[1][ipart]*momentum[1][ipart] +
                                     momentum[2][ipart]*momentum[2][ipart] );

        // Move the photons
#ifdef  __DEBUG
        for ( int i = 0 ; i<nDim_ ; i++ )
            position_old[i][ipart] = position[i][ipart];
#endif

        position[0][ipart]     += dt*momentum[0][ipart]*(*invgf)[ipart];
        position[1][ipart]     += dt*momentum[1][ipart]*(*invgf)[ipart];
        position[2][ipart]     += dt*momentum[2][ipart]*(*invgf)[ipart];
    }
}
