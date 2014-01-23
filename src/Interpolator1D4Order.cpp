
#include "Interpolator1D4Order.h" 

#include <iostream>
#include <cmath>
using namespace std;

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particle.h" 
#include "SmileiMPI_Cart1D.h" 

Interpolator1D4Order::Interpolator1D4Order(PicParams *params, SmileiMPI* smpi) : Interpolator1D(params, smpi)  
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

	dx_inv_ = 1.0/params->cell_length[0];
    
    //double defined for use in coefficients
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_48 = 1.0/48.0;
    dble_1_ov_16 = 1.0/16.0;
    dble_1_ov_12 = 1.0/12.0;
    dble_1_ov_24 = 1.0/24.0;
    dble_19_ov_96 = 19.0/96.0;
    dble_11_ov_24 = 11.0/24.0;
    dble_1_ov_384 = 1.0/384.0;
    dble_1_ov_4 = 1.0/4.0;
    dble_1_ov_6 = 1.0/6.0;
    dble_115_ov_192 = 115.0/192.0;
    dble_5_ov_8 = 5.0/8.0;

	index_domain_begin = smpi1D->getCellStartingGlobalIndex(0);

}

Interpolator1D4Order::~Interpolator1D4Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D4Order::operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc)
{
    
    // Variable declaration
    int im2, im1, i, ip1, ip2;
    double xjn, xjmxi, xjmxi2, xjmxi3, xjmxi4;
    double cim2, cim1, ci, cip1, cip2;
    
    // Static cast of the electromagnetic fields
	Field1D* Ex1D     = static_cast<Field1D*>(champs->Ex_);
	Field1D* Ey1D     = static_cast<Field1D*>(champs->Ey_);
	Field1D* Ez1D     = static_cast<Field1D*>(champs->Ez_);
	Field1D* Bx1D_m   = static_cast<Field1D*>(champs->Bx_m);
	Field1D* By1D_m   = static_cast<Field1D*>(champs->By_m);
	Field1D* Bz1D_m   = static_cast<Field1D*>(champs->Bz_m);
	
    
	// Particle position (in units of the spatial-step)
	xjn    = part->position(0)*dx_inv_;
    
    
    // --------------------------------------------------------
	// Interpolate the fields from the Primal grid : Ey, Ez, Bx
    // --------------------------------------------------------
	i      = round(xjn);      // index of the central point
	xjmxi  = xjn -(double)i;  // normalized distance to the central node
	xjmxi2 = xjmxi*xjmxi;     // square of the normalized distance to the central node
    xjmxi3 = xjmxi2*xjmxi;    // cube of the normalized distance to the central node
    xjmxi4 = xjmxi3*xjmxi;    // 4th power of the normalized distance to the central node

	// coefficients for the 4th order interpolation on 5 nodes
	cim2 = dble_1_ov_384   - dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 - dble_1_ov_12 * xjmxi3 + dble_1_ov_12 * xjmxi4;
	cim1 = dble_19_ov_96   - dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  + dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
	ci   = dble_115_ov_192 - dble_5_ov_8   * xjmxi2 + dble_1_ov_4 * xjmxi4;
    cip1 = dble_19_ov_96   + dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  - dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    cip2 = dble_1_ov_384   + dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 + dble_1_ov_12 * xjmxi3 + dble_1_ov_12 * xjmxi4;
    
	i -= index_domain_begin;
    im2    = i-2;
    im1    = i-1;
    ip1    = i+1;
    ip2    = i+2;
    
	(*ELoc).y = cim2*(*Ey1D)(im2)   + cim1*(*Ey1D)(im1)   + ci*(*Ey1D)(i)   + cip1*(*Ey1D)(ip1)   + cip2*(*Ey1D)(ip2);
	(*ELoc).z = cim2*(*Ez1D)(im2)   + cim1*(*Ez1D)(im1)   + ci*(*Ez1D)(i)   + cip1*(*Ez1D)(ip1)   + cip2*(*Ez1D)(ip2);
	(*BLoc).x = cim2*(*Bx1D_m)(im2) + cim1*(*Bx1D_m)(im1) + ci*(*Bx1D_m)(i) + cip1*(*Bx1D_m)(ip1) + cip2*(*Bx1D_m)(ip2);
	
    
    // --------------------------------------------------------
	// Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
	i      = round(xjn+0.5);  // index of the central point
	xjmxi  = xjn -(double)i;  // normalized distance to the central node
	xjmxi2 = xjmxi*xjmxi;     // square of the normalized distance to the central node
    xjmxi3 = xjmxi2*xjmxi;    // cube of the normalized distance to the central node
    xjmxi4 = xjmxi3*xjmxi;    // 4th power of the normalized distance to the central node
    
	// coefficients for the 4th order interpolation on 5 nodes
	cim2 = dble_1_ov_384   - dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 - dble_1_ov_12 * xjmxi3 + dble_1_ov_12 * xjmxi4;
	cim1 = dble_19_ov_96   - dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  + dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
	ci   = dble_115_ov_192 - dble_5_ov_8   * xjmxi2 + dble_1_ov_4 * xjmxi4;
    cip1 = dble_19_ov_96   + dble_11_ov_24 * xjmxi  + dble_1_ov_4 * xjmxi2  - dble_1_ov_6  * xjmxi3 - dble_1_ov_6  * xjmxi4;
    cip2 = dble_1_ov_384   + dble_1_ov_48  * xjmxi  + dble_1_ov_16 * xjmxi2 + dble_1_ov_12 * xjmxi3 + dble_1_ov_12 * xjmxi4;


	i -= index_domain_begin;
    im2    = i-2;
    im1    = i-1;
    ip1    = i+1;
    ip2    = i+2;

	(*ELoc).x = cim2*(*Ex1D)(im2)   + cim1*(*Ex1D)(im1)   + ci*(*Ex1D)(i)   + cip1*(*Ex1D)(ip1)   + cip2*(*Ex1D)(ip2);
	(*BLoc).y = cim2*(*By1D_m)(im2) + cim1*(*By1D_m)(im1) + ci*(*By1D_m)(i) + cip1*(*By1D_m)(ip1) + cip2*(*By1D_m)(ip2);
	(*BLoc).z = cim2*(*Bz1D_m)(im2) + cim1*(*Bz1D_m)(im1) + ci*(*Bz1D_m)(i) + cip1*(*Bz1D_m)(ip1) + cip2*(*Bz1D_m)(ip2);
	
}//END Interpolator1D4Order
