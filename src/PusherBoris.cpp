
#include "PusherBoris.h"
#include "Particle.h"

#include <iostream>
#include <cmath>

using namespace std;

PusherBoris::PusherBoris(PicParams *params, int ispec)
	: Pusher(params, ispec)
{
}

/***********************************************************************
	Lorentz Force -- leap-frog (Boris) scheme
***********************************************************************/
void PusherBoris::operator() (Particle* part, LocalFields Epart, LocalFields Bpart, double& gf)
{
	// Declaration of local variables
	// ------------------------------
	double umx, umy, umz, upx, upy, upz;
	double alpha, inv_det_T, Tx, Ty, Tz, Tx2, Ty2, Tz2;
	double TxTy, TyTz, TzTx;
	double pxsm, pysm, pzsm;
    	
	//DEBUG(5, "\tPush particle"<< part->position(0) );
	
	// --------------------------------------
	// SOLVE THE PARTICLE EQUATION OF MOTIONS
	// --------------------------------------
	
	// Half-acceleration in the electric field
	umx = part->momentum(0) + charge_over_mass_*Epart.x*dts2;
	umy = part->momentum(1) + charge_over_mass_*Epart.y*dts2;
	umz = part->momentum(2) + charge_over_mass_*Epart.z*dts2;
	gf  = sqrt( 1.0 + umx*umx + umy*umy + umz*umz );
	
	// Rotation in the magnetic field
	alpha = charge_over_mass_*dts2/gf;
	Tx    = alpha * Bpart.x;
	Ty    = alpha * Bpart.y;
	Tz    = alpha * Bpart.z;
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
	pxsm = upx + charge_over_mass_*Epart.x*dts2;
	pysm = upy + charge_over_mass_*Epart.y*dts2;
	pzsm = upz + charge_over_mass_*Epart.z*dts2;
	gf = sqrt( 1.0 + pxsm*pxsm + pysm*pysm + pzsm*pzsm );
	
	part->momentum(0) = pxsm;
	part->momentum(1) = pysm;
	part->momentum(2) = pzsm;
	
	// Move the particle
    //!\todo Make a loop on all spatial dimensions (also separate change_momentum & change_position) (MG & JD)
	part->position_old(0)  = part->position(0);
	part->position(0)     += dt*part->momentum(0)/gf;
	//DEBUG(5, "\t END "<< part->position(0) );

}

