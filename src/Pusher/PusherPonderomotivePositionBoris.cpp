#include "PusherPonderomotivePositionBoris.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;
// Pushes only position of particles interacting with envelope, not their momentum
PusherPonderomotivePositionBoris::PusherPonderomotivePositionBoris( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherPonderomotivePositionBoris::~PusherPonderomotivePositionBoris()
{
}

/**************************************************************************
    Lorentz Force + Ponderomotive force -- leap-frog (Boris-style) scheme, position advance
**************************************************************************/

void PusherPonderomotivePositionBoris::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_buffer_offset )
{

    std::vector<double> *Phi_mpart     = &( smpi->dynamics_PHI_mpart[ithread] );
    std::vector<double> *GradPhi_mpart = &( smpi->dynamics_GradPHI_mpart[ithread] );
    double *invgf = &( smpi->dynamics_invgf[ithread][0] );
    
    
    double charge_sq_over_mass_dts4, charge_sq_over_mass_sq;
    double gamma0, gamma0_sq, gamma_ponderomotive;
    double pxsm, pysm, pzsm;
    
    double* momentum_x = particles.getPtrMomentum(0);
    double* momentum_y = particles.getPtrMomentum(1);
    double* momentum_z = particles.getPtrMomentum(2);
    
    double* position_x = particles.getPtrPosition(0);
    double* position_y = NULL;
    double* position_z = NULL;
    if (nDim_>1) {
        position_y = particles.getPtrPosition(1);
        if (nDim_>2) {
            position_z = particles.getPtrPosition(2);
        }
    }
    
    short *charge = particles.getPtrCharge( ) ;
    
    int nparts;
    if (vecto) {
        nparts = GradPhi_mpart->size()/3;
    } else {
        //nparts = particles.size();
        nparts = particles.last_index.back();
    }
    
    double *Phi_m      = &( ( *Phi_mpart )[0*nparts] );
    double *GradPhi_mx = &( ( *GradPhi_mpart )[0*nparts] );
    double *GradPhi_my = &( ( *GradPhi_mpart )[1*nparts] );
    double *GradPhi_mz = &( ( *GradPhi_mpart )[2*nparts] );
    
    #ifndef _GPU
        #pragma omp simd
    #else
        int np = iend-istart;
        #pragma acc parallel present(Phi_m[istart:np], GradPhi_mx[istart:np],GradPhi_my[istart:np],GradPhi_mz[istart:np],invgf[0:nparts]) deviceptr(position_x,position_y,position_z,momentum_x,momentum_y,momentum_z,charge)
        #pragma acc loop gang worker vector
    #endif
    for( int ipart=istart ; ipart<iend; ipart++ ) { // begin loop on particles
    
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = ( double )( charge[ipart] )*( double )( charge[ipart] )*one_over_mass_*dts4;
        // (charge over mass)^2
        charge_sq_over_mass_sq      = ( double )( charge[ipart] )*one_over_mass_*( charge[ipart] )*one_over_mass_;
        
        // compute initial ponderomotive gamma
        gamma0_sq = 1.0 + momentum_x[ipart]*momentum_x[ipart] + momentum_y[ipart]*momentum_y[ipart] + momentum_z[ipart]*momentum_z[ipart] + ( *( Phi_m+ipart-ipart_buffer_offset ) )*charge_sq_over_mass_sq ;
        gamma0    = sqrt( gamma0_sq ) ;
        // ponderomotive force for ponderomotive gamma advance (Grad Phi is interpolated in time, hence the division by 2)
        pxsm = charge_sq_over_mass_dts4 * ( *( GradPhi_mx+ipart-ipart_buffer_offset ) ) / gamma0_sq ;
        pysm = charge_sq_over_mass_dts4 * ( *( GradPhi_my+ipart-ipart_buffer_offset ) ) / gamma0_sq ;
        pzsm = charge_sq_over_mass_dts4 * ( *( GradPhi_mz+ipart-ipart_buffer_offset ) ) / gamma0_sq ;
        
        // update of gamma ponderomotive
        gamma_ponderomotive = gamma0 + ( pxsm*momentum_x[ipart]+pysm*momentum_y[ipart]+pzsm*momentum_z[ipart] ) ;
        invgf[ipart-ipart_buffer_offset] = 1.0 / gamma_ponderomotive;
        
        // Move the particle
        position_x[ipart] += dt*momentum_x[ipart]*invgf[ipart-ipart_buffer_offset];
        if (nDim_>1) {
            position_y[ipart] += dt*momentum_y[ipart]*invgf[ipart-ipart_buffer_offset];
            if (nDim_>2) {
                position_z[ipart] += dt*momentum_z[ipart]*invgf[ipart-ipart_buffer_offset];
            }
        }
        
    } // end loop on particles
}
