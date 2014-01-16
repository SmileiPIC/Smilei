#include "ElectroMagn1D.h"
#include "PicParams.h"
#include "Field1D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart1D.h"

#include <iostream>
#include <math.h>
#include <sstream>
#include <string>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(PicParams* params, SmileiMPI* smpi)
	: ElectroMagn(params, smpi)
{
	// local dt to store
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	int process_coord_x = smpi1D->getProcCoord(0);
    
    oversize_ = params->oversize[0];

	// spatial-step and ratios time-step by spatial-step & spatial-step by time-step
	dx       = params->cell_length[0];
	dt       = params->timestep;
	dt_ov_dx = params->timestep/params->cell_length[0];
	dx_ov_dt = 1.0/dt_ov_dx;

	// Parameters for the Silver-Mueller boundary conditions
	Alpha_SM = 2./(1.+dt_ov_dx);
	Beta_SM  = (dt_ov_dx-1.)/(1.+dt_ov_dx);
	Gamma_SM = 4./(1.+dt_ov_dx);

	// Electromagnetic fields
	// ----------------------
	// number of nodes of the primal-grid
	nx_p = params->n_space[0]+1 + 2*params->oversize[0];
	// number of nodes of the dual-grid
	nx_d = params->n_space[0]+2 + 2*params->oversize[0];

	dimPrim.resize( params->nDim_field );
	dimDual.resize( params->nDim_field );
	
	for (size_t i=0 ; i<params->nDim_field ; i++) {
		// Standard scheme
		dimPrim[i] = params->n_space[i]+1;
		dimDual[i] = params->n_space[i]+2;
		// + Ghost domain
		dimPrim[i] += 2*params->oversize[i];
		dimDual[i] += 2*params->oversize[i];
	}
    MESSAGE( "dimPrim[0]  " <<  dimPrim[0] );

    // Allocation of the EM fields
	Ex_  = new Field1D(dimPrim, 0, false, "Ex");
	Ey_  = new Field1D(dimPrim, 1, false, "Ey");
	Ez_  = new Field1D(dimPrim, 2, false, "Ez");
	Bx_  = new Field1D(dimPrim, 0, true,  "Bx");
	By_  = new Field1D(dimPrim, 1, true,  "By");
	Bz_  = new Field1D(dimPrim, 2, true,  "Bz");
	Bx_m = new Field1D(dimPrim, 0, true,  "Bx_m");
	By_m = new Field1D(dimPrim, 1, true,  "By_m");
	Bz_m = new Field1D(dimPrim, 2, true,  "Bz_m");
	
	// Total charge currents and densities
	Jx_   = new Field1D(dimPrim, 0, false, "Jx");
	Jy_   = new Field1D(dimPrim, 1, false, "Jy");
	Jz_   = new Field1D(dimPrim, 2, false, "Jz");
	rho_  = new Field1D(dimPrim, "Rho" );
	rho_o = new Field1D(dimPrim, "Rho_old" );
	
    // Charge currents currents and density for each species
    ostringstream file_name("");
    for (unsigned int ispec=0; ispec<n_species; ispec++){
 	file_name.str("");
	file_name << "Jx_s" << ispec;
        Jx_s[ispec]  = new Field1D(dimPrim, 0, false, file_name.str().c_str());
 	file_name.str("");
	file_name << "Jy_s" << ispec;
        Jy_s[ispec]  = new Field1D(dimPrim, 1, false, file_name.str().c_str());
 	file_name.str("");
	file_name << "Jz_s" << ispec;
        Jz_s[ispec]  = new Field1D(dimPrim, 2, false, file_name.str().c_str());
 	file_name.str("");
	file_name << "rho_s" << ispec;
        rho_s[ispec] = new Field1D(dimPrim, file_name.str().c_str());
    }
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
	index_bc_min.resize( params->nDim_field, 0 );
	index_bc_max.resize( params->nDim_field, 0 );
	for (size_t i=0 ; i<params->nDim_field ; i++) {
		index_bc_min[i] = params->oversize[i];
		index_bc_max[i] = dimDual[i]-params->oversize[i]-1;
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
    MESSAGE("Entering Poisson Solver");
    
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);
	//int process_coord_x = smpi1D->getProcCoord(0);
    
    unsigned int iteration_max  = 100000;
    double       error_max      = 1.e-14; //!\todo Check what should be used to relate to machine precision
    
    double       dx_sq          = dx*dx;
    unsigned int nx_p_global    = smpi1D->n_space_global[0] + 1;
    unsigned int smilei_sz      = smpi1D->smilei_sz;
    unsigned int smilei_rk      = smpi1D->smilei_rk;
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // ----------------------------------------------------------------------------------
    vector<unsigned int> index_min_p(1);
    vector<unsigned int> index_min_d(1);
    vector<unsigned int> index_max_p(1);
    vector<unsigned int> index_max_d(1);
    index_min_p[0] = oversize_;
    index_min_d[0] = oversize_;
    index_max_p[0] = nx_p - 2 - oversize_;
    index_max_d[0] = nx_d - 2 - oversize_;
    if (smpi1D->isWester()) {
        index_min_p[0] = 0;
        index_min_d[0] = 0;
    }
    if (smpi1D->isEaster()) {
        index_max_p[0] = nx_p-1;
        index_max_d[0] = nx_d-1;
    }
    
	Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
	Field1D* rho1D = static_cast<Field1D*>(rho_);
    
    double pW=0.0;
    double pE=0.0;
    
    Field1D phi(dimPrim);    // scalar potential
    Field1D r(dimPrim);      // residual vector
    Field1D p(dimPrim);      // direction vector
    Field1D Ap(dimPrim);     // A*p vector
    
    
    // -------------------------------
    // Initialization of the variables
    // -------------------------------
    DEBUG(1,"Initialize variables");
    
    unsigned int iteration = 0;
    
    // phi: scalar potential, r: residual and p: direction
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        phi(i)   = 0.0;
        r(i)     = -dx_sq * (*rho1D)(i);
        p(i)     = r(i);
    }
    
    // scalar product of the residual (first local then reduction over all procs)
    //!\todo Generalize parallel computation of scalar product as a method for FieldND
    double rnew_dot_rnew       = 0.0;
    double rnew_dot_rnew_local = 0.0;
    //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) rnew_dot_rnew_local += r(i)*r(i);
    for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) rnew_dot_rnew_local += r(i)*r(i);
    MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    // control parameter
    //!\todo Check control parameter for // case
    double ctrl = rnew_dot_rnew/(double)(nx_p_global);
    
    
    // ---------------------------------------------------------
    // Starting iterative loop for the conjugate gradient method
    // ---------------------------------------------------------
    DEBUG(1,"Starting iterative loop");
    while ( (ctrl > error_max) && (iteration<iteration_max) ) {
        
        iteration++;
        DEBUG(5,"iteration " << iteration << " started with control parameter ctrl = " << ctrl*1.e14 << " x 1e-14");
        
        double r_dot_r = rnew_dot_rnew;
        
        // vector product Ap = A*p
        for (unsigned int i=1 ; i<dimPrim[0]-1 ; i++) Ap(i) = p(i-1) - 2.0*p(i) + p(i+1);
        
        // apply BC on Ap
        if (smpi1D->isWester()) pW        - 2.0*p(0)      + p(1);
        if (smpi1D->isEaster()) p(nx_p-2) - 2.0*p(nx_p-1) + pE;
        smpi1D->exchangeField(&Ap);
        
        // scalar product p.Ap
        double p_dot_Ap = 0.0;
        double p_dot_Ap_local = 0.0;
        //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) p_dot_Ap_local += p(i)*Ap(i);
        for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) p_dot_Ap_local += p(i)*Ap(i);
        MPI_Allreduce(&p_dot_Ap_local, &p_dot_Ap, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new potential and residual
        double alpha_k = r_dot_r/p_dot_Ap;
        for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
            phi(i) += alpha_k * p(i);
            r(i)   -= alpha_k * Ap(i);
        }
        
        // compute new residual norm
        rnew_dot_rnew       = 0.0;
        rnew_dot_rnew_local = 0.0;
        //for (unsigned int i=index_bc_min[0] ; i<index_bc_max[0]-1 ; i++) rnew_dot_rnew_local += r(i)*r(i);
        for (unsigned int i=index_min_p[0] ; i<=index_max_p[0] ; i++) rnew_dot_rnew_local += r(i)*r(i);
        MPI_Allreduce(&rnew_dot_rnew_local, &rnew_dot_rnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        // compute new direction
        double beta_k = rnew_dot_rnew/r_dot_r;
      	for (unsigned int i=0 ; i<dimPrim[0] ; i++)  p(i) = r(i) + beta_k * p(i);
        
        // compute control parameter
        ctrl = rnew_dot_rnew/(double)(nx_p_global);
        
    }//End of the iterative loop
    
    
    // Status of the solver convergence
    if (iteration == iteration_max){
        WARNING("Poisson solver did not converge: reached maximum iteration number: " << iteration
                << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    else{
        MESSAGE("Poisson solver converged at iteration: " << iteration
                << ", relative error is ctrl = " << 1.0e14*ctrl << " x 1e-14");
    }
    
    // ----------------------------------
    // Compute the electrostatic field Ex
    // ----------------------------------
    
    for (unsigned int i=1; i<nx_d-1; i++) (*Ex1D)(i) = (phi(i-1)-phi(i))/dx;
    
    
    // BC on Ex
    smpi1D->exchangeField(Ex1D);
    if (smpi1D->isWester()) (*Ex1D)(0)      = (*Ex1D)(1)      - dx*(*rho1D)(0);
    if (smpi1D->isEaster()) (*Ex1D)(nx_d-1) = (*Ex1D)(nx_d-2) + dx*(*rho1D)(nx_p-1);
    
    
    // Find field to be added to ensure BC: Ex_West = -Ex_East
    
    double Ex_West = 0.0;
    double Ex_East = 0.0;
    
    unsigned int rankWest = smpi1D->extrem_ranks[0][0];
    if (smpi1D->isWester()) {
        if (smilei_rk != smpi1D->extrem_ranks[0][0]) ERROR("western process not well defined");
        Ex_West = (*Ex1D)(index_bc_min[0]);
    }
    MPI_Bcast(&Ex_West, 1, MPI_DOUBLE, rankWest, MPI_COMM_WORLD);
    
    unsigned int rankEast = smpi1D->extrem_ranks[0][1];
    if (smpi1D->isEaster()) {
        if (smilei_rk != smpi1D->extrem_ranks[0][1]) ERROR("eastern process not well defined");
        Ex_East = (*Ex1D)(index_bc_max[0]);
    }
    MPI_Bcast(&Ex_East, 1, MPI_DOUBLE, rankEast, MPI_COMM_WORLD);
    
    
    
    cerr << "Ex_West = " << Ex_West << "  -  " << "Ex_East = " << Ex_East << endl;
    double Ex_Add = -0.5*(Ex_West+Ex_East);
    
    // Center the electrostatic field
    for (unsigned int i=0; i<nx_d; i++) (*Ex1D)(i) += Ex_Add;
    
    
    // Compute error on the Poisson equation
    double deltaPoisson_max = 0.0;
    int i_deltaPoisson_max  = -1;
    for (unsigned int i=0; i<nx_p; i++) {
        double deltaPoisson = abs( ((*Ex1D)(i+1)-(*Ex1D)(i))/dx - (*rho1D)(i) );
        if (deltaPoisson > deltaPoisson_max){
            deltaPoisson_max   = deltaPoisson;
            i_deltaPoisson_max = i;
        }
    }
    
    MESSAGE(1,"Poisson equation solved. Maximum error = " << deltaPoisson_max << " at i= " << i_deltaPoisson_max);
	
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwell(double time_dual, SmileiMPI* smpi)
{
    
    saveMagneticFields();
    solveMaxwellAmpere();
    solveMaxwellFaraday();
    applyEMBoundaryConditions(time_dual, smpi);
    smpi->exchangeB( this );
    centerMagneticFields();

}



// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::saveMagneticFields()
{
    // Static cast of the fields
	Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
	Field1D* By1D   = static_cast<Field1D*>(By_);
	Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
	Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
	Field1D* By1D_m = static_cast<Field1D*>(By_m);
	Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
    
    // for Bx^(p)
	for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
		(*Bx1D_m)(i)=(*Bx1D)(i);
	}
	//for By^(d) & Bz^(d)
	for (unsigned int i=0 ; i<dimDual[0] ; i++) {
		(*By1D_m)(i) = (*By1D)(i);
		(*Bz1D_m)(i) = (*Bz1D)(i);
	}
    
}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwellAmpere()
{

	Field1D* Ex1D = static_cast<Field1D*>(Ex_);
	Field1D* Ey1D = static_cast<Field1D*>(Ey_);
	Field1D* Ez1D = static_cast<Field1D*>(Ez_);
	Field1D* By1D = static_cast<Field1D*>(By_);
	Field1D* Bz1D = static_cast<Field1D*>(Bz_);
	Field1D* Jx1D = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D = static_cast<Field1D*>(Jz_);

	// --------------------
	// Solve Maxwell-Ampere
	// --------------------
	// Calculate the electrostatic field ex on the dual grid
	//for (unsigned int ix=0 ; ix<nx_d ; ix++){
	for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
	  (*Ex1D)(ix)= (*Ex1D)(ix) - dt* (*Jx1D)(ix) ;
	}
	// Transverse fields ey, ez  are defined on the primal grid
	//for (unsigned int ix=0 ; ix<nx_p ; ix++) {
	for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - dt * (*Jy1D)(ix) ;
		(*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - dt * (*Jz1D)(ix) ;
	}

	//smpi->exchangeE( this ); 			// Useless by construction

}


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::solveMaxwellFaraday()
{
	Field1D* Ey1D   = static_cast<Field1D*>(Ey_);
	Field1D* Ez1D   = static_cast<Field1D*>(Ez_);
	Field1D* By1D   = static_cast<Field1D*>(By_);
	Field1D* Bz1D   = static_cast<Field1D*>(Bz_);

	// ---------------------
	// Solve Maxwell-Faraday
	// ---------------------
	// NB: bx is given in 1d and defined when initializing the fields (here put to 0)
	// Transverse fields  by & bz are defined on the dual grid
	//for (unsigned int ix=1 ; ix<nx_p ; ix++) {
	for (unsigned int ix=1 ; ix<dimDual[0]-1 ; ix++) {
		(*By1D)(ix)= (*By1D)(ix) + dt_ov_dx * ( (*Ez1D)(ix) - (*Ez1D)(ix-1)) ;
		(*Bz1D)(ix)= (*Bz1D)(ix) - dt_ov_dx * ( (*Ey1D)(ix) - (*Ey1D)(ix-1)) ;
	}
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::applyEMBoundaryConditions(double time_dual, SmileiMPI* smpi)
{
	SmileiMPI_Cart1D* smpi1D = static_cast<SmileiMPI_Cart1D*>(smpi);

	Field1D* Ex1D   = static_cast<Field1D*>(Ex_);
	Field1D* Ey1D   = static_cast<Field1D*>(Ey_);
	Field1D* Ez1D   = static_cast<Field1D*>(Ez_);
	Field1D* By1D   = static_cast<Field1D*>(By_);
	Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
    
	// --------------------------------------------------
	// Laser temporal profile
	// --------------------------------------------------
	double byL=0, bzL=0, byR=0, bzR=0;

	for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
		// testing the time-profile
		// ------------------------

//		if (laser_[ilaser]->laser_struct.time_profile == "constant" ||
//            laser_[ilaser]->laser_struct.time_profile == "sin2") {
            
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

        
//		} else {
//			ERROR("Laser profile "<< ilaser <<" not allowed");
//		}//ENDif time_profile
	}//ilaser

	// ----------------------------
	// Apply EM boundary conditions
	// ----------------------------
    
	//!\todo Take care that there is a difference between primal and dual grid when putting the fields to 0 on the ghost cells (MG to JD)
    if ( smpi1D->isWester() ) {
        
		// Silver-Mueller boundary conditions (left)
		(*By1D)(0) =  Alpha_SM*(*Ez1D)(0) + Beta_SM*(*By1D)(1) + Gamma_SM*byL;
		(*Bz1D)(0) = -Alpha_SM*(*Ey1D)(0) + Beta_SM*(*Bz1D)(1) + Gamma_SM*bzL;
/*
		// Silver-Mueller boundary conditions (left)
		(*By1D)(index_bc_min[0])= Alpha_SM*(*Ez1D)(index_bc_min[0]) + Beta_SM*(*By1D)(index_bc_min[0]+1) + Gamma_SM*byL;
		(*Bz1D)(index_bc_min[0])=-Alpha_SM*(*Ey1D)(index_bc_min[0]) + Beta_SM*(*Bz1D)(index_bc_min[0]+1) + Gamma_SM*bzL;
		// Correction on unused extreme ghost
		for (unsigned int ix=0 ; ix<index_bc_min[0] ; ix++) {
			(*By1D)(ix)=0;
			(*Bz1D)(ix)=0;
			(*Ey1D)(ix)=0;
			(*Ez1D)(ix)=0;
			(*Ex1D)(ix)=0;
		}
 */
	}//if Western
    
	if ( smpi1D->isEaster() ) {
        
		// Silver-Mueller boundary conditions (right)
		(*By1D)(nx_d-1) = -Alpha_SM*(*Ez1D)(nx_p-1) + Beta_SM*(*By1D)(nx_d-2) + Gamma_SM*byR;
		(*Bz1D)(nx_d-1) =  Alpha_SM*(*Ey1D)(nx_p-1) + Beta_SM*(*Bz1D)(nx_d-2) + Gamma_SM*bzR;

/*
		// Silver-Mueller boundary conditions (right)
		(*By1D)(index_bc_max[0])=-Alpha_SM*(*Ez1D)(index_bc_max[0]) + Beta_SM*(*By1D)(index_bc_max[0]-1) + Gamma_SM*byR;
		(*Bz1D)(index_bc_max[0])= Alpha_SM*(*Ey1D)(index_bc_max[0]) + Beta_SM*(*Bz1D)(index_bc_max[0]-1) + Gamma_SM*bzR;
		// Correction on unused extreme ghost
		for (unsigned int ix=index_bc_max[0]+1 ; ix<dimDual[0] ; ix++) {
			(*By1D)(ix  )=0;
			(*Bz1D)(ix  )=0;
			(*Ex1D)(ix  )=0;
		}
		for (unsigned int ix=index_bc_max[0] ; ix<dimPrim[0] ; ix++) {
			(*Ey1D)(ix)=0;
			(*Ez1D)(ix)=0;
		}
 */
	}//if Eastern

}



// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::centerMagneticFields()
{
    // Static cast of the fields
	Field1D* Bx1D   = static_cast<Field1D*>(Bx_);
	Field1D* By1D   = static_cast<Field1D*>(By_);
	Field1D* Bz1D   = static_cast<Field1D*>(Bz_);
	Field1D* Bx1D_m = static_cast<Field1D*>(Bx_m);
	Field1D* By1D_m = static_cast<Field1D*>(By_m);
	Field1D* Bz1D_m = static_cast<Field1D*>(Bz_m);
    
    // for Bx^(p)
	for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
		(*Bx1D_m)(i) = ( (*Bx1D)(i)+ (*Bx1D_m)(i))*0.5 ;
    }
    
	// for By^(d) & Bz^(d)
	for (unsigned int i=0 ; i<dimDual[0] ; i++) {
		(*By1D_m)(i)= ((*By1D)(i)+(*By1D_m)(i))*0.5 ;
		(*Bz1D_m)(i)= ((*Bz1D)(i)+(*Bz1D_m)(i))*0.5 ;
    }
    
}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::restartRhoJ()
{
	Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);
	Field1D* rho1D_o = static_cast<Field1D*>(rho_o);
	
    // --------------------------
    // Total currents and density
    // --------------------------
    
	// put longitudinal current to zero on the dual grid
    for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
		(*Jx1D)(ix)    = 0.0;
	}
    
	// all fields are defined on the primal grid
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*rho1D_o)(ix) = (*rho1D)(ix);
		(*rho1D)(ix)   = 0.0;
		(*Jy1D)(ix)    = 0.0;
		(*Jz1D)(ix)    = 0.0;
	}
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++){
        Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
        
        // put longitudinal current to zero on the dual grid
        for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
            (*Jx1D_s)(ix)  = 0.0;
        }
        
        // all fields are defined on the primal grid
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
            (*rho1D_s)(ix) = 0.0;
            (*Jy1D_s)(ix)  = 0.0;
            (*Jz1D_s)(ix)  = 0.0;
        }
    }//END loop on species ispec
}



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalRhoJ()
{
	Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
	Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
	Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
	Field1D* rho1D   = static_cast<Field1D*>(rho_);

    for (unsigned int ispec=0; ispec<n_species; ispec++){
        Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
        
        // put longitudinal current to zero on the dual grid
        for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
            (*Jx1D)(ix)  += (*Jx1D_s)(ix);
        }
        
        // all fields are defined on the primal grid
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
            (*Jy1D)(ix)  += (*Jy1D_s)(ix);
            (*Jz1D)(ix)  += (*Jz1D_s)(ix);
            (*rho1D)(ix) += (*rho1D_s)(ix);
        }
    }//END loop on species ispec
}

/*
// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge density and transverse currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::initRho()
{
    MESSAGE(" >>> DIMPRIM >>> DIMDUAL >>>" << dimPrim[0] << " " << dimDual[0]);
	Field1D* rho1D_o    = static_cast<Field1D*>(rho1D_o);
    Field1D* rho1D      = static_cast<Field1D*>(rho1D);
    
	// all fields are defined on the primal grid
	//for (unsigned int i=0 ; i<nx_p ; i++) {
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
		(*rho1D_o)(ix) = 0.0;
		(*rho1D)(ix)   = 0.0;
	}
    
}
*/
