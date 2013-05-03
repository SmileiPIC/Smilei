
#include "ElectroMagn1D.h"

#include "PicParams.h"
#include "Field1D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"

#include <iostream>
#include <math.h>
using namespace std;

ElectroMagn1D::ElectroMagn1D(PicParams* params, SmileiMPI* smpi)
  : ElectroMagn(params, smpi)
{
	dimPrim.resize( params->nDim_field );
	dimDual.resize( params->nDim_field );	
	dspace.resize   ( params->nDim_field );
	dspacesdt.resize( params->nDim_field );
	dtsdspace.resize( params->nDim_field );
	
	for (size_t i=0 ; i<params->nDim_field ; i++) {
		dimPrim[i] = params->n_space[i]+1;
		dimDual[i] = params->n_space[i]+2;
		dspace[i] = params->cell_length[i];
		dspacesdt[i] = params->cell_length[i]/params->timestep;
		dtsdspace[i] = params->timestep/params->cell_length[i];
	}
	// Parameters for the Silver-Mueller boundary conditions 
	A_ = 4./(1.+dtsdspace[0]);
	B_ = (dtsdspace[0]-1.)/(1.+dtsdspace[0]);
	C_ = 2./(1.+dtsdspace[0]);
	
	// Laser carac To Do
	// a0_delta_y_L_;
	// a0_delta_z_L_;
	// a0_delta_y_R_;
	// a0_delta_z_R_;
	// tau1_L_, tau2_L_, tau1_R_, tau2_R_;
	// las_tordr_L_, las_tordr_R_;
	
	Ex_ = new Field1D( dimDual, "fex" );
	Ey_ = new Field1D( dimPrim, "fey" );
	Ez_ = new Field1D( dimPrim, "fez" );
	
	Bx_ = new Field1D( dimPrim, "fbx" );
	By_ = new Field1D( dimDual, "fby" );
	Bz_ = new Field1D( dimDual, "fbz" );
	
	Bx_m = new Field1D(dimPrim);
	By_m = new Field1D(dimDual);
	Bz_m = new Field1D(dimDual);
	
	Jx_ = new Field1D(dimDual, "fjx");
	Jy_ = new Field1D(dimPrim, "fjy");
	Jz_ = new Field1D(dimPrim, "fjz");
	
	rho_ = new Field1D(dimPrim, "rho");
	rho_o = new Field1D(dimPrim, "rho_old");

	iPrim_beg.resize(params->nDim_field, 0);
	iPrim_end.resize(params->nDim_field, 0);
	iDual_beg.resize(params->nDim_field, 0);
	iDual_end.resize(params->nDim_field, 0);
	iDual_beg_nobc.resize(params->nDim_field, 0);
	iDual_end_nobc.resize(params->nDim_field, 0);

	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	for (size_t i=0 ; i<params->nDim_field ; i++) {
	  iPrim_beg[0] = params->oversize[0];
	  iPrim_end[0] = dimPrim[0] - params->oversize[0];
	  iDual_beg[0] = params->oversize[0];
	  iDual_end[0] = dimDual[0] - params->oversize[0];
	  if ( smpi1D->isWester() ) iDual_beg_nobc[0] = iDual_beg[0]+1;
	  else                      iDual_beg_nobc[0] = iDual_beg[0];
	  if ( smpi1D->isEaster() ) iDual_end_nobc[0] = iDual_end[0]-1;
	  else                      iDual_end_nobc[0] = iDual_end[0];

	}
		
}//END constructor Electromagn1D


ElectroMagn1D::~ElectroMagn1D()
{
}


void ElectroMagn1D::initMaxwell(SmileiMPI* smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
	Field1D* rho1D = static_cast<Field1D*>(rho_);

	if (process_coord_x==0) (*Ex1D)(iDual_beg[0]) = 0.0;
	for ( unsigned int ix = iDual_beg_nobc[0] ; ix < iDual_end_nobc[0] ; ix++ ){
	  (*Ex1D)(ix) = (*Ex1D)(ix-1) + dspace[0]* (*rho1D)(ix-1);
	}

}//END initMaxwell


/***********************************************************************************************************************
 Explicit (leap-frog) Solver of the Maxwell-Faraday-Ampere 
 Equations are discretized on a Yee Mesh -- FDTD scheme
 **********************************************************************************************************************/
void ElectroMagn1D::solveMaxwell(double time_dual, double dt, SmileiMPI* smpi)
{
	Field1D* Ex1D = static_cast<Field1D*>(Ex_);
	Field1D* Ey1D = static_cast<Field1D*>(Ey_);
	Field1D* Ez1D = static_cast<Field1D*>(Ez_);
	Field1D* Bx1D = static_cast<Field1D*>(Bx_);
	Field1D* By1D = static_cast<Field1D*>(By_);
	Field1D* Bz1D = static_cast<Field1D*>(Bz_);
	Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
	Field1D* By1D_m = static_cast<Field1D*>(By_m);
	Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
	Field1D* Jx1D = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D = static_cast<Field1D*>(Jz_);

	DEBUG(5,"solveMaxwell ElectroMagn1D " << time_dual);

	// Define laser-polarization
	double byL=0, bzL=0, byR=0, bzR=0;
	
	for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {		
		if (laser_[ilaser]->laser_struct.time_profile == "constant") {
			if (laser_[ilaser]->laser_struct.angle == 0) {
				// Incident field (left boundary)
				byL += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
				bzL += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
			} else if (laser_[ilaser]->laser_struct.angle == 180) {
				// Incident field (right boundary)
				byR += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);				
				bzR += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
            } else {
				ERROR("Angle not allowed for 1D laser pulse " << ilaser);
			}
		} else {
			ERROR("Laser profile "<< ilaser <<" not allowed");
		}
	}


	// SAVE FORMER B-FIELD TO CENTER THE B-FIELD
	// IN THE PARTICLE PUSHER
	// AND OF THE E-FIELD TO CENTER FOR DIAG
	// -----------------------------------------
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*Bx1D_m)(ix)=(*Bx1D)(ix);
	}
	for (unsigned int ix=0 ; ix<dimDual[0] ; ix++){
		(*By1D_m)(ix)=(*By1D)(ix);
		(*Bz1D_m)(ix)=(*Bz1D)(ix);
	}

	
	// SOLVE MAXWELL-AMPERE
	//---------------------
	// Calculate the electrostatic field ex
	for (unsigned int ix=iDual_beg[0] ; ix<iDual_end[0] ; ix++) {
		(*Ex1D)(ix)= (*Ex1D)(ix) - dt* (*Jx1D)(ix) ;
	}
	
	// Transverse fields ey, ez  are defined on the primal grid
	for (unsigned int ix=iPrim_beg[0] ; ix<iPrim_end[0] ; ix++) {
		(*Ey1D)(ix)= (*Ey1D)(ix) - dtsdspace[0] * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - dt * (*Jy1D)(ix) ;
		(*Ez1D)(ix)= (*Ez1D)(ix) + dtsdspace[0] * ( (*By1D)(ix+1) - (*By1D)(ix)) - dt * (*Jz1D)(ix) ;
	}
	smpi->exchangeE( this );


	// SOLVE MAXWELL-FARADAY
	// ---------------------
	// NB: bx is given in 1d ==> defined in init_fields (here put to 0)  
	// Transverse fields  by & bz are defined on the dual grid
	for (unsigned int ix=iDual_beg_nobc[0] ; ix<iDual_end_nobc[0] ; ix++) {
		(*By1D)(ix)= (*By1D)(ix) + dtsdspace[0] * ( (*Ez1D)(ix) - (*Ez1D)(ix-1)) ;
		(*Bz1D)(ix)= (*Bz1D)(ix) - dtsdspace[0] * ( (*Ey1D)(ix) - (*Ey1D)(ix-1)) ;
	}
	smpi->exchangeB( this );

	// BOUNDARY CONDITIONS
	// -------------------
	if ( iDual_beg_nobc[0] != iDual_beg[0] ) {
	  // Silver-Mueller boundary conditions (left)
	  (*By1D)(iDual_beg[0])= A_*byL + B_* (*By1D)(iDual_beg[0]+1) + C_* (*Ez1D)(iDual_beg[0]); 
	  (*Bz1D)(iDual_beg[0])= A_*bzL + B_* (*Bz1D)(iDual_beg[0]+1) - C_* (*Ey1D)(iDual_beg[0]);
	}
	if ( iDual_end_nobc[0] != iDual_end[0] ) {
	  // Silver-Mueller boundary conditions (right)
	  (*By1D)(iDual_end[0])= A_*byR + B_* (*By1D)(iDual_end[0]-1) - C_* (*Ez1D)(iDual_end[0]);
	  (*Bz1D)(iDual_end[0])= A_*bzR + B_* (*Bz1D)(iDual_end[0]-1) + C_* (*Ey1D)(iDual_end[0]);
	}
	
	// CENTER THE B-FIELD FOR THE PARTICLE PUSHER
	// ------------------------------------------
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++)
		(*Bx1D_m)(ix)= ( (*Bx1D)(ix)+ (*Bx1D_m)(ix))*0.5 ;
	for (unsigned int ix=0 ; ix<dimDual[0] ; ix++){
		(*By1D_m)(ix)= ((*By1D)(ix)+(*By1D_m)(ix))*0.5 ;
		(*Bz1D_m)(ix)= ((*Bz1D)(ix)+(*Bz1D_m)(ix))*0.5 ;
	}
	
}


/***********************************************************************************************************************
 Calculate the longitudinal current by solving the charge-conservation equation
 **********************************************************************************************************************/
void ElectroMagn1D::chargeConserving(SmileiMPI* smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	Field1D* Jx1D = static_cast<Field1D*>(Jx_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);

	if (process_coord_x==0) (*Jx1D)(iDual_beg[0]) = 0.0;
	for ( unsigned int ix = iDual_beg_nobc[0] ; ix < iDual_end_nobc[0] ; ix++ ){
	  (*Jx1D)(ix) = (*Jx1D)(ix-1) - dspacesdt[0] * ( (*rho1D)(ix-1)-(*rho1D_o)(ix-1) );
	}

}


/***********************************************************************************************************************
 Reinitialize the total charge density and transverse currents
 - save current density as old density (charge conserving scheme)
 - put the new density and currents to 0
 **********************************************************************************************************************/
void ElectroMagn1D::initRhoJ()
{
	Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);
	
	//defined on the primal grid
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*rho1D_o)(ix) = (*rho1D)(ix);
		(*rho1D)(ix)   = 0.0;
		(*Jy1D)(ix)    = 0.0;
		(*Jz1D)(ix)    = 0.0;
	}
    
}
