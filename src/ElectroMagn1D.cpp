#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Field1D.h"
#include "Laser.h"

#include <iostream>
#include <math.h>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(PicParams* params) : ElectroMagn(params)
{
    // number of nodes of the primal-grid
    nx_p = params->n_space[0]+1;
    
    // number of nodes of the dual-grid
    nx_d = params->n_space[0]+2;
    
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
    std::vector<unsigned int> dimPrim; dimPrim.resize(1); dimPrim[0] = params->n_space[0]+1;
    std::vector<unsigned int> dimDual; dimDual.resize(1); dimDual[0] = params->n_space[0]+2;
    
	Ex_  = new Field1D( dimDual, "fex" );
    Ey_  = new Field1D( dimPrim, "fey" );
	Ez_  = new Field1D( dimPrim, "fez" );
	Bx_  = new Field1D( dimPrim, "fbx" );
	By_  = new Field1D( dimDual, "fby" );
	Bz_  = new Field1D( dimDual, "fbz" );
	Bx_m = new Field1D( dimPrim );
	By_m = new Field1D( dimDual );
	Bz_m = new Field1D( dimDual );
	
    // Total charge currents and densities
	Jx_   = new Field1D( dimDual, "fjx");
	Jy_   = new Field1D( dimPrim, "fjy");
	Jz_   = new Field1D( dimPrim, "fjz");
	rho_  = new Field1D( dimPrim, "rho");
	rho_o = new Field1D( dimPrim, "rho_old");
	
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
void ElectroMagn1D::solvePoisson()
{
	Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
	Field1D* rho1D = static_cast<Field1D*>(rho_);
	
	// Initialize the electrostatic field by solving Poisson at t = 0
    // \todo Generalise this so one minimises the electrostatic energy (MG)
	(*Ex1D)(0) = 0.0;
	for (unsigned int i=1 ; i<nx_d ; i++)
    {
		(*Ex1D)(i) = (*Ex1D)(i-1) + dx * (*rho1D)(i-1);
    }
	
}//END solvePoisson



// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwell(double time_dual, double dt)
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
	Field1D* Jx1D   = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D   = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D   = static_cast<Field1D*>(Jz_);
	
	DEBUG(5,"solveMaxwell ElectroMagn1D " << time_dual);
	
    
    // --------------------------------------------------
	// Define the laser fields at left & right boundaries
    // --------------------------------------------------
    
	double byL=0, bzL=0, byR=0, bzR=0;
	
	for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++)
    {
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
	for (unsigned int i=0 ; i<nx_p ; i++){
		(*Bx1D_m)(i) = (*Bx1D)(i);
	}
	for (unsigned int i=0 ; i<nx_d ; i++){
		(*By1D_m)(i) = (*By1D)(i);
		(*Bz1D_m)(i) = (*Bz1D)(i);
	}
	
	
    // --------------------
	// Solve Maxwell-Ampere
	// --------------------
	
	// Calculate the electrostatic field ex on the dual grid
	for (unsigned int i=0 ; i<nx_d ; i++){
		(*Ex1D)(i) = (*Ex1D)(i) - dt * (*Jx1D)(i) ;
    }
	
	// Transverse fields ey, ez  are defined on the primal grid
	for (unsigned int i=0 ; i<nx_p ; i++) {
		(*Ey1D)(i) = (*Ey1D)(i) - dt_ov_dx * ( (*Bz1D)(i+1) - (*Bz1D)(i)) - dt * (*Jy1D)(i) ;
		(*Ez1D)(i) = (*Ez1D)(i) + dt_ov_dx * ( (*By1D)(i+1) - (*By1D)(i)) - dt * (*Jz1D)(i) ;
	}
	
    
    // ---------------------
	// Solve Maxwell-Faraday
	// ---------------------
    
	// NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    
	// Transverse fields  by & bz are defined on the dual grid
	for (unsigned int i=1 ; i<nx_p ; i++) {
		(*By1D)(i) = (*By1D)(i) + dt_ov_dx * ( (*Ez1D)(i) - (*Ez1D)(i-1)) ;
		(*Bz1D)(i) = (*Bz1D)(i) - dt_ov_dx * ( (*Ey1D)(i) - (*Ey1D)(i-1)) ;
	}
	
    
    // ----------------------------
	// Apply EM boundary conditions
	// ----------------------------
    
    //!\todo Make boundary conditions on the EM fields as an external method (MG)
    
	// Silver-Mueller boundary conditions (left)
	(*By1D)(0) = A_*byL + B_* (*By1D)(1) + C_* (*Ez1D)(0) ;
	(*Bz1D)(0) = A_*bzL + B_* (*Bz1D)(1) - C_* (*Ey1D)(0) ;
	
    // Silver-Mueller boundary conditions (right)
	(*By1D)(nx_d-1) = A_*byR + B_* (*By1D)(nx_d-2) - C_* (*Ez1D)(nx_p-1) ;
	(*Bz1D)(nx_d-1) = A_*bzR + B_* (*Bz1D)(nx_d-2) + C_* (*Ey1D)(nx_p-1) ;
	
    
    // ------------------------------------------------
	// Center the magnetic fields (for particle pusher)
	// ------------------------------------------------
    
	for (unsigned int i=0 ; i<nx_p ; i++){
		(*Bx1D_m)(i) = ( (*Bx1D)(i)+ (*Bx1D_m)(i))*0.5;
    }
	for (unsigned int i=0 ; i<nx_d ; i++){
		(*By1D_m)(i) = ((*By1D)(i)+(*By1D_m)(i))*0.5;
		(*Bz1D_m)(i) = ((*Bz1D)(i)+(*Bz1D_m)(i))*0.5;
	}
	
}


// ---------------------------------------------------------------------------------------------------------------------
// Calculate the longitudinal current by solving the charge-conservation equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::chargeConserving()
{
	Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);
	
    //!\todo Replace this by Esirkepov method for the calculation of longitudinal currents (MG)
    // longitudinal currents defined on the dual-grid
	(*Jx1D)(0)=0.0;
	for (unsigned int i=1 ; i<nx_d ; i++) {
		(*Jx1D)(i) = (*Jx1D)(i-1) - dx_ov_dt * ( (*rho1D)(i-1)-(*rho1D_o)(i-1) );
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
    for (unsigned int i=0 ; i<nx_d ; i++){
		(*Jx1D)(i)    = 0.0;
	}
    
	// all fields are defined on the primal grid
	for (unsigned int i=0 ; i<nx_p ; i++)
    {
		(*rho1D_o)(i) = (*rho1D)(i);
		(*rho1D)(i)   = 0.0;
		(*Jy1D)(i)    = 0.0;
		(*Jz1D)(i)    = 0.0;
	}
    
}
