#include "Interpolator1D3Order.h"

#include <cmath>
#include <iostream>

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particles.h"
#include "SmileiMPI_Cart1D.h"

using namespace std;

Interpolator1D3Order::Interpolator1D3Order(PicParams &params, SmileiMPI *smpi) : Interpolator1D(params, smpi)
{
    dx_inv_   = 1.0/params.cell_length[0];

    dble_1ov6 = 1.0/6.0;
    dble_2ov3 = 2.0/3.0;
}

/***********************************************************************
	Interpolate the field fx defined on the primal grid
	with size nstp_x and space step stp_x_inv at the position
	xj and return the value fxj
***********************************************************************/
void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc)
{
    double xjn, xi, xi2, xi3;

    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D* Ex1D     = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(EMfields->Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(EMfields->By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(EMfields->Bz_m);


    // Calculate the normalized positions
    // ----------------------------------
    xjn    = particles.position(0, ipart)*dx_inv_;


    // Primal Grid : Ey, Ez, Bx
    // ------------------------
    ip_ = (int)xjn;            // index of the 2nd node
    xi  = xjn - (double)ip_;   //normalized distance to the 2nd node
    xi2 = xi*xi;
    xi3 = xi*xi*xi;

    // 3rd order interpolation on 4 nodes
    coeffp_[0]  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffp_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffp_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffp_[3]  = xi3*dble_1ov6;

    ip_ -= index_domain_begin;

    (*ELoc).y = compute(coeffp_, Ey1D,   ip_);  
    (*ELoc).z = compute(coeffp_, Ez1D,   ip_);  
    (*BLoc).x = compute(coeffp_, Bx1D_m, ip_);

    // Dual Grid : Ex, By, Bz
    // ----------------------
    id_   = (int)(xjn+0.50);        // position of the 2nd node
    xi  = xjn - (double)id_ + 0.5;    // normalized distance to the central node
    xi2 = xi*xi;
    xi3 = xi*xi3;

    // 3rd order interpolation on 4 nodes
    coeffd_[0] = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
    coeffd_[1]  = dble_2ov3 - xi2 + 0.5*xi3;
    coeffd_[2]  = dble_1ov6 + 0.5*(xi+xi2-xi3);
    coeffd_[3]  = xi3*dble_1ov6;

    id_ -= index_domain_begin;

    (*ELoc).x = compute(coeffd_, Ex1D,   id_);  
    (*BLoc).y = compute(coeffd_, By1D_m, id_);  
    (*BLoc).z = compute(coeffd_, Bz1D_m, id_);  

}//END Interpolator1D3Order


void Interpolator1D3Order::operator() (ElectroMagn* EMfields, Particles &particles, int ipart, LocalFields* ELoc, LocalFields* BLoc, LocalFields* JLoc, double* RhoLoc){

    // Interpolate E, B
    // Compute coefficient for ipart position    (*this)(EMfields, particles, ipart, ELoc, BLoc);

    //!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
    Field1D* Jx1D     = static_cast<Field1D*>(EMfields->Jx_);
    Field1D* Jy1D     = static_cast<Field1D*>(EMfields->Jy_);
    Field1D* Jz1D     = static_cast<Field1D*>(EMfields->Jz_);
    Field1D* Rho1D    = static_cast<Field1D*>(EMfields->rho_);
    
    // Primal Grid : Jy, Jz, Rho
    // ------------------------
    (*JLoc).y = compute(coeffp_, Jy1D,  ip_);  
    (*JLoc).z = compute(coeffp_, Jz1D,  ip_);  
    (*RhoLoc) = compute(coeffp_, Rho1D, ip_);    
    
    // Dual Grid : Jx
    // ----------------------
    (*JLoc).x = compute(coeffd_, Jx1D,  id_);  
    
}
