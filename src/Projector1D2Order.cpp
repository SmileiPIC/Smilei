
#include "Projector1D2Order.h" 

#include <iostream>
#include <cmath>
using namespace std;

#include "ElectroMagn.h"
#include "Field1D.h"
#include "Particle.h" 

#include "Tools.h"

/***********************************************************************************************************************
 Projection of the total density and currents on the primal grid
 **********************************************************************************************************************/
Projector1D2Order::Projector1D2Order (PicParams* params) :Projector1D(params)
{
	dx_inv_ = 1.0/params->cell_length[0];               // inverse of the spatial-step
	DEBUG("cell_length "<< params->cell_length[0]);
}


void Projector1D2Order::operator() (ElectroMagn* champs, Particle* part, double gf)
{
	Field1D* rho1D  = static_cast<Field1D*>(champs->rho_);
	Field1D* Jy1D   = static_cast<Field1D*>(champs->Jy_);	
	Field1D* Jz1D   = static_cast<Field1D*>(champs->Jz_);

	
    //Declaration of local variables
    int i;
	double xjn,xjmxi,xjmxi2;
    double rho_j = part->chargeDensity();  // charge density of the macro-particle
    double cry_j = rho_j*part->moments(1)/gf;       // current density allow the y-direction of the macroparticle
    double crz_j = rho_j*part->moments(2)/gf;       // current density allow the y-direction of the macroparticle
    

	//Locate particle on the primal grid
	xjn    = part->position(0) * dx_inv_;     // normalized distance to the first node
	i      = round(xjn);                  // index of the central node
	xjmxi  = xjn - (double)i;             // normalized distance to the nearest grid point
	xjmxi2 = xjmxi*xjmxi;                 // square of the normalized distance to the nearest grid point

	
    // 2nd order projection for the total density
	(*rho1D)( i-1) = ((*rho1D)(i-1) + 0.5 * (xjmxi2-xjmxi+0.25) * rho_j );
	(*rho1D)( i  ) = ((*rho1D)(i  ) +  (0.75-xjmxi2)            * rho_j );
	(*rho1D)( i+1) = ((*rho1D)(i+1) + 0.5 * (xjmxi2+xjmxi+0.25) * rho_j );

    // 2nd order projection for the total current
	(*Jy1D)( i-1) = ((*Jy1D)(i-1) + 0.5 * (xjmxi2-xjmxi+0.25) * cry_j );
	(*Jy1D)( i  ) = ((*Jy1D)(i  ) +  (0.75-xjmxi2)            * cry_j );
	(*Jy1D)( i+1) = ((*Jy1D)(i+1) + 0.5 * (xjmxi2+xjmxi+0.25) * cry_j );

	(*Jz1D)( i-1) = ((*Jz1D)(i-1) + 0.5 * (xjmxi2-xjmxi+0.25) * crz_j );
	(*Jz1D)( i  ) = ((*Jz1D)(i  ) +  (0.75-xjmxi2)            * crz_j );
	(*Jz1D)( i+1) = ((*Jz1D)(i+1) + 0.5 * (xjmxi2+xjmxi+0.25) * crz_j );
	
}

void Projector1D2Order::operator() (Field* rho, Particle* part)
{
	Field1D* rho1D  = static_cast<Field1D*>(rho);

	
	//Declaration of local variables
	int i;
	double xjn,xjmxi,xjmxi2;
	double rho_j = part->chargeDensity();  // charge density of the macro-particle
    

	//Locate particle on the grid
	xjn    = part->position(0) * dx_inv_;     // normalized distance to the first node
	i      = round(xjn);                  // index of the central node
	xjmxi  = xjn - (double)i;             // normalized distance to the nearest grid point
	xjmxi2 = xjmxi*xjmxi;                 // square of the normalized distance to the nearest grid point

	
	// 2nd order projection for the total density
	(*rho1D)( i-1)  = ((*rho1D)(i-1) + 0.5 * (xjmxi2-xjmxi+0.25) * rho_j );
	(*rho1D)( i  )  = ((*rho1D)(i  ) +  (0.75-xjmxi2)            * rho_j );
	(*rho1D)( i+1)  = ((*rho1D)(i+1) + 0.5 * (xjmxi2+xjmxi+0.25) * rho_j );

}

