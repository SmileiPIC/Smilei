
#include "Interpolator1D3Order.h" 

#include <iostream>
#include <cmath>
using namespace std;

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particle.h" 
#include "SmileiMPI_Cart1D.h" 

Interpolator1D3Order::Interpolator1D3Order(PicParams *params, SmileiMPI* smpi) : Interpolator1D(params, smpi)  
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

	dx_inv_   = 1.0/params->cell_length[0];
	dble_1ov6 = 1.0/6.0;
	dble_2ov3 = 2.0/3.0;

	int process_coord_x = smpi1D->getProcCoord(0);
	index_domain_begin =  process_coord_x*(params->n_space[0]-2*params->oversize[0]);
}

/***********************************************************************
	Interpolate the field fx defined on the primal grid 
	with size nstp_x and space step stp_x_inv at the position 
	xj and return the value fxj
***********************************************************************/
void Interpolator1D3Order::operator() (ElectroMagn* champs, Particle* part, chLocaux* ELoc, chLocaux* BLoc)
{
	int i, im1, ip1, ip2;
	double xjn, xi, xi2, xi3;
	double c1, c2, c3, c4;
    
	//!\todo Julien, can you check that this is indeed the centered B-field which is passed to the pusher?
	Field1D* Ex1D     = static_cast<Field1D*>(champs->Ex_);
	Field1D* Ey1D     = static_cast<Field1D*>(champs->Ey_);
	Field1D* Ez1D     = static_cast<Field1D*>(champs->Ez_);
	Field1D* Bx1D_m   = static_cast<Field1D*>(champs->Bx_m);
	Field1D* By1D_m   = static_cast<Field1D*>(champs->By_m);
	Field1D* Bz1D_m   = static_cast<Field1D*>(champs->Bz_m);
	
    
	// Calculate the normalized positions
	// ----------------------------------
	xjn    = part->position(0)*dx_inv_;
    
    
	// Primal Grid : Ey, Ez, Bx
	// ------------------------
	i = (int)xjn;            // index of the 2nd node
	im1 = i-1;
	ip1 = i+1;
	ip2 = i+2;
	xi  = xjn - (double)i;   //normalized distance to the 2nd node
	xi2 = xi*xi;
	xi3 = xi*xi*xi;
    
	// 3rd order interpolation on 4 nodes
	c1  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
	c2  = dble_2ov3 - xi2 + 0.5*xi3;
	c3  = dble_1ov6 + 0.5*(xi+xi2-xi3);
	c4  = xi3*dble_1ov6;

	i -= index_domain_begin;
    
	(*ELoc).y = c1 * (*Ey1D)(im1)   + c2 * (*Ey1D)(i)   + c3 * (*Ey1D)(ip1)   + c4 * (*Ey1D)(ip2);
	(*ELoc).z = c1 * (*Ez1D)(im1)   + c2 * (*Ez1D)(i)   + c3 * (*Ez1D)(ip1)   + c4 * (*Ez1D)(ip2);
	(*BLoc).x = c1 * (*Bx1D_m)(im1) + c2 * (*Bx1D_m)(i) + c3 * (*Bx1D_m)(ip1) + c4 * (*Bx1D_m)(ip2);
	
    
	// Dual Grid : Ex, By, Bz
	// ----------------------
	i   = (int)(xjn+0.50);        // position of the 2nd node
	im1 = i-1;
	ip1 = i+1;
	ip2 = i+2;
	xi  = xjn - (double)i + 0.5;   // normalized distance to the central node
	xi2 = xi*xi;
	xi3 = xi*xi3;
	
	// 3rd order interpolation on 4 nodes
	c1  = (1.0-xi3)*dble_1ov6 - 0.5*(xi-xi2);
	c2  = dble_2ov3 - xi2 + 0.5*xi3;
	c3  = dble_1ov6 + 0.5*(xi+xi2-xi3);
	c4  = xi3*dble_1ov6;

	i -= index_domain_begin;
	
	(*ELoc).x = c1 * (*Ex1D)(im1)   + c2 * (*Ex1D)(i)   + c3 * (*Ex1D)(ip1)   + c4 * (*Ex1D)(ip2);
	(*BLoc).y = c1 * (*By1D_m)(im1) + c2 * (*By1D_m)(i) + c3 * (*By1D_m)(ip1) + c4 * (*By1D_m)(ip2);
	(*BLoc).z = c1 * (*Bz1D_m)(im1) + c2 * (*Bz1D_m)(i) + c3 * (*Bz1D_m)(ip1) + c4 * (*Bz1D_m)(ip2);
	
}//END Interpolator1D3Order

