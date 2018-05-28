#include "PusherPonderomotivePositionBoris.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;
// Pushes only position of particles interacting with envelope, not their momentum
PusherPonderomotivePositionBoris::PusherPonderomotivePositionBoris(Params& params, Species *species)
    : Pusher(params, species)
{
}

PusherPonderomotivePositionBoris::~PusherPonderomotivePositionBoris()
{
}

/**************************************************************************
    Lorentz Force + Ponderomotive force -- leap-frog (Boris-style) scheme, position advance
**************************************************************************/

void PusherPonderomotivePositionBoris::operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ipart_ref)
{
    std::vector<double> *Phipart        = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPhipart    = &(smpi->dynamics_GradPHIpart[ithread]);
    std::vector<double> *Phioldpart     = &(smpi->dynamics_PHIoldpart[ithread]);
    std::vector<double> *GradPhioldpart = &(smpi->dynamics_GradPHIoldpart[ithread]);
    
    double charge_sq_over_mass_dts4,charge_over_mass_sq;
    double inv_gamma0,inv_gamma_ponderomotive;
    double pxsm, pysm, pzsm;
    
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
    
    short* charge = &( particles.charge(0) );
    
    int nparts = particles.size();
    
    double* Phi         = &( (*Phipart)[0*nparts] );
    double* Phiold      = &( (*Phioldpart)[0*nparts] );
    double* GradPhix    = &( (*GradPhipart)[0*nparts] );
    double* GradPhiy    = &( (*GradPhipart)[1*nparts] );
    double* GradPhiz    = &( (*GradPhipart)[2*nparts] );
    double* GradPhioldx = &( (*GradPhioldpart)[0*nparts] );
    double* GradPhioldy = &( (*GradPhioldpart)[1*nparts] );
    double* GradPhioldz = &( (*GradPhioldpart)[2*nparts] );
    
    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) { // begin loop on particles
    
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = (double)(charge[ipart])*(double)(charge[ipart])*one_over_mass_*dts4;         
        // (charge over mass)^2
        charge_over_mass_sq      = (double)(charge[ipart])*one_over_mass_*(charge[ipart])*one_over_mass_;

        // compute initial ponderomotive gamma (more precisely, its inverse) 
        inv_gamma0 = 1./sqrt( 1. + momentum[0][ipart]*momentum[0][ipart] + momentum[1][ipart]*momentum[1][ipart] + momentum[2][ipart]*momentum[2][ipart] + (*(Phi+ipart)+*(Phiold+ipart))*charge_over_mass_sq*0.5 );
    
        // ponderomotive force for ponderomotive gamma advance (Grad Phi is interpolated in time, hence the division by 2)
        pxsm = charge_sq_over_mass_dts4 * ( *(GradPhix+ipart) + *(GradPhioldx+ipart) ) * 0.5 * inv_gamma0 ;
        pysm = charge_sq_over_mass_dts4 * ( *(GradPhiy+ipart) + *(GradPhioldy+ipart) ) * 0.5 * inv_gamma0 ;
        pzsm = charge_sq_over_mass_dts4 * ( *(GradPhiz+ipart) + *(GradPhioldz+ipart) ) * 0.5 * inv_gamma0 ;
    
        // update of gamma ponderomotive (more precisely, the inverse)
        inv_gamma_ponderomotive = 1./( 1./inv_gamma0 + (pxsm*momentum[0][ipart]+pysm*momentum[1][ipart]+pzsm*momentum[2][ipart])*inv_gamma0 );
  
        // Move the particle
#ifdef  __DEBUG
        for ( int i = 0 ; i<nDim_ ; i++ ) 
          position_old[i][ipart] = position[i][ipart];
#endif
        for ( int i = 0 ; i<nDim_ ; i++ ) 
            position[i][ipart]     += dt*momentum[i][ipart]*inv_gamma_ponderomotive;

    } // end loop on particles
}
