#include "PusherPonderomotiveBorisV.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;
// Pushes only momentum of particles interacting with envelope, not their position
PusherPonderomotiveBorisV::PusherPonderomotiveBorisV(Params& params, Species *species)
    : Pusher(params, species)
{
}

PusherPonderomotiveBorisV::~PusherPonderomotiveBorisV()
{
}

/**************************************************************************
    Lorentz Force + Ponderomotive force -- leap-frog (Boris-style) scheme
**************************************************************************/

void PusherPonderomotiveBorisV::operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)
{

    std::vector<double> *Epart       = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart       = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *Phipart     = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPhipart = &(smpi->dynamics_GradPHIpart[ithread]);

    double charge_over_mass_dts2,charge_sq_over_mass_dts4;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double local_invgf;
    double inv_gamma0,inv_gamma_ponderomotive;
    double charge_sq_over_mass_sq;

    int* cell_keys;

    double* momentum[3];
    for ( int i = 0 ; i<3 ; i++ )
        momentum[i] =  &( particles.momentum(i,0) );

    short* charge = &( particles.charge(0) );

    int nparts = particles.size();
    double* Ex       = &( (*Epart)[0*nparts] );
    double* Ey       = &( (*Epart)[1*nparts] );
    double* Ez       = &( (*Epart)[2*nparts] );
    double* Bx       = &( (*Bpart)[0*nparts] );
    double* By       = &( (*Bpart)[1*nparts] );
    double* Bz       = &( (*Bpart)[2*nparts] );
    double* Phi      = &( (*Phipart)[0*nparts] );
    double* GradPhix = &( (*GradPhipart)[0*nparts] );
    double* GradPhiy = &( (*GradPhipart)[1*nparts] );
    double* GradPhiz = &( (*GradPhipart)[2*nparts] );
  
    particles.cell_keys.resize(nparts);
    cell_keys = &( particles.cell_keys[0]);

    double dcharge[nparts];
    #pragma omp simd
    for (int ipart=0 ; ipart<nparts; ipart++ ) 
        dcharge[ipart] = (double)(charge[ipart]);

    #pragma omp simd
    for (int ipart=0 ; ipart<nparts; ipart++ ) {
        double psm[3], um[3];

        charge_over_mass_dts2 = dcharge[ipart]*one_over_mass_*dts2;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = dcharge[ipart]*dcharge[ipart]*one_over_mass_*dts4;         
        // (charge over mass)^2
        charge_sq_over_mass_sq   = dcharge[ipart]*one_over_mass_*dcharge[ipart]*one_over_mass_;


        // compute initial ponderomotive gamma (more precisely, its inverse) 
        inv_gamma0 = 1./sqrt( 1. + momentum[0][ipart]*momentum[0][ipart] + momentum[1][ipart]*momentum[1][ipart] + momentum[2][ipart]*momentum[2][ipart] + *(Phi+ipart)*charge_sq_over_mass_sq );
        
        // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
        psm[0] = inv_gamma0 * (charge_over_mass_dts2*(*(Ex+ipart)) - charge_sq_over_mass_dts4*(*(GradPhix+ipart)) * inv_gamma0 ) * momentum[0][ipart];
        psm[1] = inv_gamma0 * (charge_over_mass_dts2*(*(Ey+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiy+ipart)) * inv_gamma0 ) * momentum[1][ipart];
        psm[2] = inv_gamma0 * (charge_over_mass_dts2*(*(Ez+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiz+ipart)) * inv_gamma0 ) * momentum[2][ipart];
        
        // update of gamma ponderomotive (more precisely, the inverse)
        inv_gamma_ponderomotive = 1./( 1./inv_gamma0 + (psm[0]+psm[1]+psm[2])*0.5 );

        // init Half-acceleration in the electric field and ponderomotive force 
        psm[0] = charge_over_mass_dts2 * (*(Ex+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhix+ipart)) * inv_gamma_ponderomotive ;
        psm[1] = charge_over_mass_dts2 * (*(Ey+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhiy+ipart)) * inv_gamma_ponderomotive ;
        psm[2] = charge_over_mass_dts2 * (*(Ez+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhiz+ipart)) * inv_gamma_ponderomotive ;

        um[0] = momentum[0][ipart] + psm[0];
        um[1] = momentum[1][ipart] + psm[1];
        um[2] = momentum[2][ipart] + psm[2];

        // Rotation in the magnetic field, using updated gamma ponderomotive
        alpha = charge_over_mass_dts2 * inv_gamma_ponderomotive;
        Tx    = alpha * (*(Bx+ipart));
        Ty    = alpha * (*(By+ipart));
        Tz    = alpha * (*(Bz+ipart));
        Tx2   = Tx*Tx;
        Ty2   = Ty*Ty;
        Tz2   = Tz*Tz;
        TxTy  = Tx*Ty;
        TyTz  = Ty*Tz;
        TzTx  = Tz*Tx;
        inv_det_T = 1.0/(1.0+Tx2+Ty2+Tz2);

        psm[0] += (  (1.0+Tx2-Ty2-Tz2)* um[0]  +      2.0*(TxTy+Tz)* um[1]  +      2.0*(TzTx-Ty)* um[2]  )*inv_det_T;
        psm[0] += (      2.0*(TxTy-Tz)* um[0]  +  (1.0-Tx2+Ty2-Tz2)* um[1]  +      2.0*(TyTz+Tx)* um[2]  )*inv_det_T;
        psm[0] += (      2.0*(TzTx+Ty)* um[0]  +      2.0*(TyTz-Tx)* um[1]  +  (1.0-Tx2-Ty2+Tz2)* um[2]  )*inv_det_T;

        momentum[0][ipart] = psm[0];
        momentum[1][ipart] = psm[1];
        momentum[2][ipart] = psm[2];

    }

}
