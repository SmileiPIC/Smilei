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

    std::vector<double> *Phi_mpart     = &(smpi->dynamics_PHI_mpart[ithread]);
    std::vector<double> *GradPhi_mpart = &(smpi->dynamics_GradPHI_mpart[ithread]);
    std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);

    
    double charge_sq_over_mass_dts4,charge_sq_over_mass_sq;
    double gamma0,gamma0_sq,gamma_ponderomotive;
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
    
    double* Phi_m      = &( (*Phi_mpart)[0*nparts] );
    double* GradPhi_mx = &( (*GradPhi_mpart)[0*nparts] );
    double* GradPhi_my = &( (*GradPhi_mpart)[1*nparts] );
    double* GradPhi_mz = &( (*GradPhi_mpart)[2*nparts] );
    
    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) { // begin loop on particles
    
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = (double)(charge[ipart])*(double)(charge[ipart])*one_over_mass_*dts4;         
        // (charge over mass)^2
        charge_sq_over_mass_sq      = (double)(charge[ipart])*one_over_mass_*(charge[ipart])*one_over_mass_;

        // compute initial ponderomotive gamma 
        gamma0_sq = 1. + momentum[0][ipart]*momentum[0][ipart] + momentum[1][ipart]*momentum[1][ipart] + momentum[2][ipart]*momentum[2][ipart] + (*(Phi_m+ipart))*charge_sq_over_mass_sq ;
        gamma0    = sqrt(gamma0_sq) ;
        // ponderomotive force for ponderomotive gamma advance (Grad Phi is interpolated in time, hence the division by 2)
        pxsm = charge_sq_over_mass_dts4 * ( *(GradPhi_mx+ipart) ) / gamma0_sq ;
        pysm = charge_sq_over_mass_dts4 * ( *(GradPhi_my+ipart) ) / gamma0_sq ;
        pzsm = charge_sq_over_mass_dts4 * ( *(GradPhi_mz+ipart) ) / gamma0_sq ;
    
        // update of gamma ponderomotive 
        gamma_ponderomotive = gamma0 + (pxsm*momentum[0][ipart]+pysm*momentum[1][ipart]+pzsm*momentum[2][ipart]) ;
        (*invgf)[ipart] = 1. / gamma_ponderomotive;
  
        // Move the particle
#ifdef  __DEBUG
        for ( int i = 0 ; i<nDim_ ; i++ ) 
          position_old[i][ipart] = position[i][ipart];
#endif
        for ( int i = 0 ; i<nDim_ ; i++ ) 
            position[i][ipart]     += dt*momentum[i][ipart]/gamma_ponderomotive;

    } // end loop on particles
}
