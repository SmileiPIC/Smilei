#include "PusherBoris.h"

#include <iostream>
#include <cmath>

#include "Species.h"

#include "Particles.h"

using namespace std;

PusherBoris::PusherBoris(Params& params, Species *species)
    : Pusher(params, species)
{
}

PusherBoris::~PusherBoris()
{
}

/***********************************************************************
	Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/

void PusherBoris::operator() (Particles &particles, SmileiMPI* smpi, int istart, int iend, int ithread)
{
    std::vector<LocalFields> *Epart = &(smpi->dynamics_Epart[ithread]);
    std::vector<LocalFields> *Bpart = &(smpi->dynamics_Bpart[ithread]);
    std::vector<double> *gf = &(smpi->dynamics_gf[ithread]);

    double charge_over_mass_ ;
    double umx, umy, umz, upx, upy, upz;
    double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
    double TxTy, TyTz, TzTx;
    double pxsm, pysm, pzsm;

    for (int ipart=istart ; ipart<iend; ipart++ ) {
        charge_over_mass_ = static_cast<double>(particles.charge(ipart))*one_over_mass_;
        //(*this)(particles, ipart, (*Epart)[ipart], (*Bpart)[ipart] , (*gf)[ipart]);
        umx = particles.momentum(0, ipart) + charge_over_mass_*(*Epart)[ipart].x*dts2;
        umy = particles.momentum(1, ipart) + charge_over_mass_*(*Epart)[ipart].y*dts2;
        umz = particles.momentum(2, ipart) + charge_over_mass_*(*Epart)[ipart].z*dts2;
        (*gf)[ipart]  = sqrt( 1.0 + umx*umx + umy*umy + umz*umz );

        // Rotation in the magnetic field
        alpha = charge_over_mass_*dts2/(*gf)[ipart];
        Tx    = alpha * (*Bpart)[ipart].x;
        Ty    = alpha * (*Bpart)[ipart].y;
        Tz    = alpha * (*Bpart)[ipart].z;
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

        // Half-acceleration in the electric field
        pxsm = upx + charge_over_mass_*(*Epart)[ipart].x*dts2;
        pysm = upy + charge_over_mass_*(*Epart)[ipart].y*dts2;
        pzsm = upz + charge_over_mass_*(*Epart)[ipart].z*dts2;
        (*gf)[ipart] = sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );

        particles.momentum(0, ipart) = pxsm;
        particles.momentum(1, ipart) = pysm;
        particles.momentum(2, ipart) = pzsm;

        // Move the particle
        for ( int i = 0 ; i<nDim_ ; i++ ) 
            particles.position(i, ipart)     += dt*particles.momentum(i, ipart)/(*gf)[ipart];
    }
}
