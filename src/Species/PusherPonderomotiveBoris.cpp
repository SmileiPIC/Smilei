#include "PusherPonderomotiveBoris.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;
// Pushes only momentum of particles interacting with envelope, not their position
PusherPonderomotiveBoris::PusherPonderomotiveBoris(Params& params, Species *species)
    : Pusher(params, species)
{
}

PusherPonderomotiveBoris::~PusherPonderomotiveBoris()
{
}

/**************************************************************************
    Lorentz Force + Ponderomotive force -- leap-frog (Boris-style) scheme, momentum advance
**************************************************************************/

void PusherPonderomotiveBoris::operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread, int ipart_ref)
{
    std::vector<double> *Epart       = &(smpi->dynamics_Epart[ithread]);
    std::vector<double> *Bpart       = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *Phipart     = &(smpi->dynamics_PHIpart[ithread]);
    std::vector<double> *GradPhipart = &(smpi->dynamics_GradPHIpart[ithread]);


    //std::vector<double> *invgf = &(smpi->dynamics_invgf[ithread]);
    double charge_over_mass_dts2,charge_sq_over_mass_dts4;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;
    double local_invgf;
    double gamma0,gamma0_sq,gamma_ponderomotive;
    double charge_sq_over_mass_sq;

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

    #pragma omp simd
    for (int ipart=istart ; ipart<iend; ipart++ ) {

        charge_over_mass_dts2    = (double)(charge[ipart])*one_over_mass_*dts2;
        // ! ponderomotive force is proportional to charge squared and the field is divided by 4 instead of 2
        charge_sq_over_mass_dts4 = (double)(charge[ipart])*(double)(charge[ipart])*one_over_mass_*dts4;         
        // (charge over mass)^2
        charge_sq_over_mass_sq   = (double)(charge[ipart])*one_over_mass_*(charge[ipart])*one_over_mass_;

        // compute initial ponderomotive gamma 
        gamma0_sq  = 1. + momentum[0][ipart]*momentum[0][ipart] + momentum[1][ipart]*momentum[1][ipart] + momentum[2][ipart]*momentum[2][ipart] + *(Phi+ipart)*charge_sq_over_mass_sq ;
        gamma0     = sqrt(gamma0_sq) ;
        // ( electric field + ponderomotive force for ponderomotive gamma advance ) scalar multiplied by momentum
        pxsm = (gamma0*charge_over_mass_dts2*(*(Ex+ipart)) - charge_sq_over_mass_dts4*(*(GradPhix+ipart))  ) * momentum[0][ipart] / gamma0_sq;
        pysm = (gamma0*charge_over_mass_dts2*(*(Ey+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiy+ipart))  ) * momentum[1][ipart] / gamma0_sq;
        pzsm = (gamma0*charge_over_mass_dts2*(*(Ez+ipart)) - charge_sq_over_mass_dts4*(*(GradPhiz+ipart))  ) * momentum[2][ipart] / gamma0_sq;
        
        // update of gamma ponderomotive 
        gamma_ponderomotive = gamma0 + (pxsm+pysm+pzsm)*0.5 ;

        // init Half-acceleration in the electric field and ponderomotive force 
        pxsm = charge_over_mass_dts2 * (*(Ex+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhix+ipart)) / gamma_ponderomotive ;
        pysm = charge_over_mass_dts2 * (*(Ey+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhiy+ipart)) / gamma_ponderomotive ;
        pzsm = charge_over_mass_dts2 * (*(Ez+ipart)) - charge_sq_over_mass_dts4 * (*(GradPhiz+ipart)) / gamma_ponderomotive ;

        umx = momentum[0][ipart] + pxsm;
        umy = momentum[1][ipart] + pysm;
        umz = momentum[2][ipart] + pzsm;

        // Rotation in the magnetic field, using updated gamma ponderomotive
        alpha = charge_over_mass_dts2 / gamma_ponderomotive;
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

        upx = (  (1.0+Tx2-Ty2-Tz2)* umx  +      2.0*(TxTy+Tz)* umy  +      2.0*(TzTx-Ty)* umz  )*inv_det_T;
        upy = (      2.0*(TxTy-Tz)* umx  +  (1.0-Tx2+Ty2-Tz2)* umy  +      2.0*(TyTz+Tx)* umz  )*inv_det_T;
        upz = (      2.0*(TzTx+Ty)* umx  +      2.0*(TyTz-Tx)* umy  +  (1.0-Tx2-Ty2+Tz2)* umz  )*inv_det_T;

        // finalize Half-acceleration in the electric field and ponderomotive force
        pxsm += upx;
        pysm += upy;
        pzsm += upz;
        
        momentum[0][ipart] = pxsm;
        momentum[1][ipart] = pysm;
        momentum[2][ipart] = pzsm;

    }
}
