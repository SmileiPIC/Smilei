#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Field1D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"

#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(PicParams* params, SmileiMPI* smpi)
	: ElectroMagn(params, smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	// spatial-step and ratios time-step by spatial-step & spatial-step by time-step
	dx       = params->cell_length[0];
	dt_ov_dx = params->timestep/params->cell_length[0];
	dx_ov_dt = 1.0/dt_ov_dx;

	// Parameters for the Silver-Mueller boundary conditions
	A_ = 4./(1.+dt_ov_dx);
	B_ = (dt_ov_dx-1.)/(1.+dt_ov_dx);
	C_ = 2./(1.+dt_ov_dx);

	// Electromagnetic fields
	// ----------------------
	// number of nodes of the primal-grid
	nx_p = params->n_space[0]+1;
	// number of nodes of the dual-grid
	nx_d = params->n_space[0]+2;

	//std::vector<unsigned int> dimPrim;
	//std::vector<unsigned int> dimDual;

	dimPrim.resize( params->nDim_field );
	dimDual.resize( params->nDim_field );	
	
	for (size_t i=0 ; i<params->nDim_field ; i++) {
	  dimPrim[i] = params->n_space[i]+1; // +1 really used on rank 0
	  dimDual[i] = params->n_space[i]+2; // +1 needed but +2 to preserve scheme
	}
	cout << process_coord_x << " : " << dimPrim[0] << " " << dimDual[0] << endl;

	ostringstream name;
	name.str(""); name << "fex." << process_coord_x;
	Ex_ = new Field1D( dimDual, name.str() );
	name.str(""); name << "fey." << process_coord_x;
	Ey_ = new Field1D( dimPrim, name.str() );
	name.str(""); name << "fez." << process_coord_x;
	Ez_ = new Field1D( dimPrim, name.str() );
	name.str(""); name << "fbx." << process_coord_x;
	Bx_ = new Field1D( dimPrim, name.str() );
	name.str(""); name << "fby." << process_coord_x;
	By_ = new Field1D( dimDual, name.str() );
	name.str(""); name << "fbz." << process_coord_x;
	Bz_ = new Field1D( dimDual, name.str() );
	Bx_m = new Field1D(dimPrim);
	By_m = new Field1D(dimDual);
	Bz_m = new Field1D(dimDual);
	
	// Total charge currents and densities
	name.str(""); name << "fjx." << process_coord_x;
	Jx_ = new Field1D(dimDual, name.str() );
	name.str(""); name << "fjy." << process_coord_x;
	Jy_ = new Field1D(dimPrim, name.str() );
	name.str(""); name << "fjz." << process_coord_x;
	Jz_ = new Field1D(dimPrim, name.str() );
	name.str(""); name << "rho." << process_coord_x;
	rho_ = new Field1D(dimPrim, name.str() );
	name.str(""); name << "rho_old." << process_coord_x;
	rho_o = new Field1D(dimPrim, name.str() );

	iPrim_beg.resize(params->nDim_field, 0);
	iPrim_end.resize(params->nDim_field, 0);
	iDual_beg.resize(params->nDim_field, 0);
	iDual_end.resize(params->nDim_field, 0);
	iDual_beg_nobc.resize(params->nDim_field, 0);
	iDual_end_nobc.resize(params->nDim_field, 0);


	for (size_t i=0 ; i<params->nDim_field ; i++) {
	  iPrim_beg[0] = params->oversize[0];
	  iPrim_end[0] = dimPrim[0] - params->oversize[0];
	  iDual_beg[0] = params->oversize[0];
	  iDual_end[0] = dimDual[0] - params->oversize[0];

	  //if ( !smpi1D->isWester() ) iDual_beg[0]++;
	  //if (  smpi1D->isWester() ) iDual_end[0]--;

	  if ( smpi1D->isWester() ) iDual_beg_nobc[0] = iDual_beg[0]+1;
	  else                      iDual_beg_nobc[0] = iDual_beg[0];
	  if ( smpi1D->isEaster() ) iDual_end_nobc[0] = iDual_end[0]-1;
	  else                      iDual_end_nobc[0] = iDual_end[0];

	}
		
}//END constructor Electromagn1D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::~ElectroMagn1D()
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solvePoisson(SmileiMPI* smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
	Field1D* rho1D = static_cast<Field1D*>(rho_);

	// Initialize the electrostatic field by solving Poisson at t = 0
	// \todo Generalise this so one minimises the electrostatic energy (MG)
        //if (process_coord_x==0) (*Ex1D)(iDual_beg[0]) = 0.0;
	//for (unsigned int i=1 ; i<nx_d ; i++) {
	for ( unsigned int ix = iDual_beg[0] ; ix < iDual_end[0] ; ix++ ){
	  //if (process_coord_x==0)  
	    (*Ex1D)(ix) = (*Ex1D)(ix-1) + dx* (*rho1D)(ix-1);
	    //if (process_coord_x!=0)  (*Ex1D)(ix) = (*Ex1D)(ix-1) + dx* (*rho1D)(ix);
	}
	
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwell(double time_dual, double dt, SmileiMPI* smpi)
{
	Field1D* Ex1D   = static_cast<Field1D*>(Ex_);
	Field1D* Ey1D   = static_cast<Field1D*>(Ey_);
	Field1D* Ez1D   = static_cast<Field1D*>(Ez_);
	Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
	Field1D* By1D   = static_cast<Field1D*>(By_);
	Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
	Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
	Field1D* By1D_m = static_cast<Field1D*>(By_m);
	Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
	Field1D* Jx1D = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D = static_cast<Field1D*>(Jz_);

	DEBUG(5,"solveMaxwell ElectroMagn1D " << time_dual);

	// --------------------------------------------------
	// Define the laser fields at left & right boundaries
	// --------------------------------------------------

	double byL=0, bzL=0, byR=0, bzR=0;
	
	for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
		// testing the time-profile
		// ------------------------
        
		if (laser_[ilaser]->laser_struct.time_profile == "constant") {
			if (laser_[ilaser]->laser_struct.angle == 0){
				// Incident field (left boundary)
				byL += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
				bzL += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
			} else if (laser_[ilaser]->laser_struct.angle == 180){
				// Incident field (right boundary)
				byR += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);				
				bzR += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
			} else {
				ERROR("Angle not allowed for 1D laser pulse " << ilaser);
			}
		} else {
			ERROR("Laser profile "<< ilaser <<" not allowed");
		}//ENDif time_profile
	}//ilaser

	// ----------------------------------------------
	// Save the magnetic fields (used to center them)
	// ----------------------------------------------
	//for (unsigned int ix=0 ; ix<nx_p ; ix++) {
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*Bx1D_m)(ix)=(*Bx1D)(ix);
	}
	//for (unsigned int ix=0 ; ix<nx_d ; ix++) {
	for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
		(*By1D_m)(ix) = (*By1D)(ix);
		(*Bz1D_m)(ix) = (*Bz1D)(ix);
	}

	// --------------------
	// Solve Maxwell-Ampere
	// --------------------
	// Calculate the electrostatic field ex on the dual grid
	//for (unsigned int ix=0 ; ix<nx_d ; ix++){
	//for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
	for (unsigned int ix=iDual_beg[0] ; ix<iDual_end[0] ; ix++) {
	  (*Ex1D)(ix)= (*Ex1D)(ix) - dt* (*Jx1D)(ix) ;
	}
	// Transverse fields ey, ez  are defined on the primal grid
	//for (unsigned int ix=0 ; ix<nx_p ; ix++) {
	for (unsigned int ix=iPrim_beg[0] ; ix<iPrim_end[0] ; ix++) {
		(*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - dt * (*Jy1D)(ix) ;
		(*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - dt * (*Jz1D)(ix) ;
	}
	smpi->exchangeE( this );

	// ---------------------
	// Solve Maxwell-Faraday
	// ---------------------
	// NB: bx is given in 1d and defined when initializing the fields (here put to 0)
	// Transverse fields  by & bz are defined on the dual grid
	//for (unsigned int ix=1 ; ix<nx_p ; ix++) {
	for (unsigned int ix=iDual_beg_nobc[0] ; ix<iDual_end_nobc[0] ; ix++) {
		(*By1D)(ix)= (*By1D)(ix) + dt_ov_dx * ( (*Ez1D)(ix) - (*Ez1D)(ix-1)) ;
		(*Bz1D)(ix)= (*Bz1D)(ix) - dt_ov_dx * ( (*Ey1D)(ix) - (*Ey1D)(ix-1)) ;
	}
	smpi->exchangeB( this );

	// ----------------------------
	// Apply EM boundary conditions
	// ----------------------------
	//!\todo Make boundary conditions on the EM fields as an external method (MG)
	if ( iDual_beg_nobc[0] != iDual_beg[0] ) {
	  // Silver-Mueller boundary conditions (left)
	  (*By1D)(iDual_beg[0])= A_*byL + B_* (*By1D)(iDual_beg[0]+1) + C_* (*Ez1D)(iDual_beg[0]); 
	  (*Bz1D)(iDual_beg[0])= A_*bzL + B_* (*Bz1D)(iDual_beg[0]+1) - C_* (*Ey1D)(iDual_beg[0]);
	}
	if ( iDual_end_nobc[0] != iDual_end[0] ) {
	  // Silver-Mueller boundary conditions (right)
	  (*By1D)(iDual_end_nobc[0])= A_*byR + B_* (*By1D)(iDual_end_nobc[0]-1) - C_* (*Ez1D)(iDual_end_nobc[0]);
	  (*Bz1D)(iDual_end_nobc[0])= A_*bzR + B_* (*Bz1D)(iDual_end_nobc[0]-1) + C_* (*Ey1D)(iDual_end_nobc[0]);
	}
	{
		/*
		// Silver-Mueller boundary conditions (left)
		(*By1D)(0) = A_*byL + B_* (*By1D)(1) + C_* (*Ez1D)(0) ;
		(*Bz1D)(0) = A_*bzL + B_* (*Bz1D)(1) - C_* (*Ey1D)(0) ;
		// Silver-Mueller boundary conditions (right)
		(*By1D)(nx_d-1) = A_*byR + B_* (*By1D)(nx_d-2) - C_* (*Ez1D)(nx_p-1) ;
		(*Bz1D)(nx_d-1) = A_*bzR + B_* (*Bz1D)(nx_d-2) + C_* (*Ey1D)(nx_p-1) ;
		*/
	}

	// ------------------------------------------------
	// Center the magnetic fields (for particle pusher)
	// ------------------------------------------------
	//for (unsigned int ix=0 ; ix<nx_p ; ix++)
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++)
		(*Bx1D_m)(ix)= ( (*Bx1D)(ix)+ (*Bx1D_m)(ix))*0.5 ;
	//for (unsigned int ix=0 ; ix<nx_d ; ix++) {
	for (unsigned int ix=0 ; ix<dimDual[0] ; ix++){
		(*By1D_m)(ix)= ((*By1D)(ix)+(*By1D_m)(ix))*0.5 ;
		(*Bz1D_m)(ix)= ((*Bz1D)(ix)+(*Bz1D_m)(ix))*0.5 ;
    	}
	
}


// ---------------------------------------------------------------------------------------------------------------------
// Calculate the longitudinal current by solving the charge-conservation equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::chargeConserving(SmileiMPI* smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);

	Field1D* Jx1D = static_cast<Field1D*>(Jx_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);

	//!\todo Replace this by Esirkepov method for the calculation of longitudinal currents (MG)
	// longitudinal currents defined on the dual-grid
	if (process_coord_x==0) (*Jx1D)(iDual_beg[0]) = 0.0;
	//for (unsigned int i=1 ; i<nx_d ; i++) {
	for ( unsigned int ix = iDual_beg_nobc[0] ; ix < iDual_end_nobc[0] ; ix++ ){
		(*Jx1D)(ix) = (*Jx1D)(ix-1) - dx_ov_dt * ( (*rho1D)(ix-1)-(*rho1D_o)(ix-1) );
	}

}


// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::initRhoJ()
{
	Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);
	
	// put longitudinal current to zero on the dual grid
	//for (unsigned int i=0 ; i<nx_d ; i++){
        for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
		(*Jx1D)(ix)    = 0.0;
	}
    
	// all fields are defined on the primal grid
	//for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*rho1D_o)(ix) = (*rho1D)(ix);
		(*rho1D)(ix)   = 0.0;
		(*Jy1D)(ix)    = 0.0;
		(*Jz1D)(ix)    = 0.0;
	}
    
}
