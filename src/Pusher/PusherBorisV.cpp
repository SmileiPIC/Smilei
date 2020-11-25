#include "PusherBorisV.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherBorisV::PusherBorisV( Params &params, Species *species )
    : Pusher( params, species )
{
}

PusherBorisV::~PusherBorisV()
{
}

/***********************************************************************
    Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBorisV::operator()( Particles &particles, SmileiMPI *smpi, int istart, int iend, int ithread, int ipart_ref )
{
    std::vector<double> *Epart = &( smpi->dynamics_Epart[ithread] );
    std::vector<double> *Bpart = &( smpi->dynamics_Bpart[ithread] );
    double *invgf = &( smpi->dynamics_invgf[ithread][0] );
    
    double charge_over_mass_dts2;
    double inv_det_T, Tx, Ty, Tz;
    double local_invgf;
    //int IX;
    
    //int* cell_keys;
    
    double *momentum[3];
    for( int i = 0 ; i<3 ; i++ ) {
        momentum[i] =  &( particles.momentum( i, 0 ) );
    }
    double *position[3];
    for( int i = 0 ; i<nDim_ ; i++ ) {
        position[i] =  &( particles.position( i, 0 ) );
    }
    short *charge = &( particles.charge( 0 ) );
    
    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    double *Bx = &( ( *Bpart )[0*nparts] );
    double *By = &( ( *Bpart )[1*nparts] );
    double *Bz = &( ( *Bpart )[2*nparts] );
    
    //particles.cell_keys.resize(nparts);
    //cell_keys = &( particles.cell_keys[0]);
    
    vector<double> dcharge(nparts);
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        dcharge[ipart-ipart_ref] = ( double )( charge[ipart] );
    }
    
    #pragma omp simd
    for( int ipart=istart ; ipart<iend; ipart++ ) {
        double psm[3], um[3];
        
        charge_over_mass_dts2 = dcharge[ipart-ipart_ref]*one_over_mass_*dts2;
        
        // init Half-acceleration in the electric field
        psm[0] = charge_over_mass_dts2*( *( Ex+ipart-ipart_ref ) );
        psm[1] = charge_over_mass_dts2*( *( Ey+ipart-ipart_ref ) );
        psm[2] = charge_over_mass_dts2*( *( Ez+ipart-ipart_ref ) );
        
        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*invgf)[ipart]);
        um[0] = momentum[0][ipart] + psm[0];
        um[1] = momentum[1][ipart] + psm[1];
        um[2] = momentum[2][ipart] + psm[2];
        
        // Rotation in the magnetic field
        local_invgf = charge_over_mass_dts2 / sqrt( 1.0 + um[0]*um[0] + um[1]*um[1] + um[2]*um[2] );
        Tx    = local_invgf * ( *( Bx+ipart-ipart_ref ) );
        Ty    = local_invgf * ( *( By+ipart-ipart_ref ) );
        Tz    = local_invgf * ( *( Bz+ipart-ipart_ref ) );
        inv_det_T = 1.0/( 1.0+Tx*Tx+Ty*Ty+Tz*Tz );
        
        psm[0] += ( ( 1.0+Tx*Tx-Ty*Ty-Tz*Tz )* um[0]  +      2.0*( Tx*Ty+Tz )* um[1]  +      2.0*( Tz*Tx-Ty )* um[2] )*inv_det_T;
        psm[1] += ( 2.0*( Tx*Ty-Tz )* um[0]  + ( 1.0-Tx*Tx+Ty*Ty-Tz*Tz )* um[1]  +      2.0*( Ty*Tz+Tx )* um[2] )*inv_det_T;
        psm[2] += ( 2.0*( Tz*Tx+Ty )* um[0]  +      2.0*( Ty*Tz-Tx )* um[1]  + ( 1.0-Tx*Tx-Ty*Ty+Tz*Tz )* um[2] )*inv_det_T;
        
        // finalize Half-acceleration in the electric field
        local_invgf = 1. / sqrt( 1.0 + psm[0]*psm[0] + psm[1]*psm[1] + psm[2]*psm[2] );
        invgf[ipart-ipart_ref] = local_invgf;
        
        momentum[0][ipart] = psm[0];
        momentum[1][ipart] = psm[1];
        momentum[2][ipart] = psm[2];
        
        // Move the particle
        local_invgf *= dt;
        for( int i = 0 ; i<nDim_ ; i++ ) {
            position[i][ipart]     += psm[i]*local_invgf;
        }
        
    }
    
    // This is temporarily moved to SpeciesV.cpp
    //#pragma omp simd
    //for (int ipart=istart ; ipart<iend; ipart++ )  {
    //
    //    for ( int i = 0 ; i<nDim_ ; i++ ){
    //        cell_keys[ipart] *= nspace[i];
    //        cell_keys[ipart] += round( (position[i][ipart]-min_loc_vec[i]) * dx_inv_[i] );
    //    }
    //
    //}
    
}
