
#include "Interpolator1D2Order.h" 

#include <iostream>
#include <cmath>
using namespace std;

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particle.h" 

// ---------------------------------------------------------------------------------------------------------------------
// Creator for Interpolator1D2Order
// ---------------------------------------------------------------------------------------------------------------------
Interpolator1D2Order::Interpolator1D2Order(PicParams *params) : Interpolator1D(params)  
{
	dx_inv_ = 1.0/params->cell_length[0];
}



// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator1D2Order::operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc)
{
    
    // Variable declaration
    int i;
    double xjn, xjmxi, xjmxi2;
    double coeffInf, coeffCur, coeffSup;
    
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
	xjmxi2 = pow(xjmxi,2);    // square of the normalized distance to the central node
	
    // 2nd order interpolation on 3 nodes
	coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
	coeffCur = (0.75-xjmxi2);
	coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);
    
	(*ELoc).y =  coeffInf * (*Ey1D)(i-1)   + coeffCur * (*Ey1D)(i)   + coeffSup * (*Ey1D)(i+1);
	(*ELoc).z =  coeffInf * (*Ez1D)(i-1)   + coeffCur * (*Ez1D)(i)   + coeffSup * (*Ez1D)(i+1);
	(*BLoc).x =  coeffInf * (*Bx1D_m)(i-1) + coeffCur * (*Bx1D_m)(i) + coeffSup * (*Bx1D_m)(i+1);
	
    
    // --------------------------------------------------------
	// Interpolate the fields from the Dual grid : Ex, By, Bz
    // --------------------------------------------------------
	i      = round(xjn+0.5);        // index of the central point
	xjmxi  = xjn - (double)i +0.5;  // normalized distance to the central node
	xjmxi2 = pow(xjmxi,2);          // square of the normalized distance to the central node
	
    // 2nd order interpolation on 3 nodes
	coeffInf = 0.5 * (xjmxi2-xjmxi+0.25);
	coeffCur = (0.75-xjmxi2);
	coeffSup = 0.5 * (xjmxi2+xjmxi+0.25);
	
	(*ELoc).x =  coeffInf * (*Ex1D)(i-1)   + coeffCur * (*Ex1D)(i)   + coeffSup * (*Ex1D)(i+1);
	(*BLoc).y =  coeffInf * (*By1D_m)(i-1) + coeffCur * (*By1D_m)(i) + coeffSup * (*By1D_m)(i+1);
	(*BLoc).z =  coeffInf * (*Bz1D_m)(i-1) + coeffCur * (*Bz1D_m)(i) + coeffSup * (*Bz1D_m)(i+1);
	
}//END Interpolator1D2Order
