#include "ElectroMagn2D.h"
#include "PicParams.h"
#include "Field2D.h"
#include "Laser.h"

#include "SmileiMPI.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <math.h>
#include <sstream>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(PicParams* params, SmileiMPI* smpi)
    : ElectroMagn(params, smpi)
{
    // local dt to store
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
	int process_coord_x = smpi2D->getProcCoord(0);
    
    
    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------
    
    // time-step
    dt       = params->timestep;
    
	// spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
	dx       = params->cell_length[0];
	dt_ov_dx = dt/dx;
	dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
	dy       = params->cell_length[1];
	dt_ov_dy = dt/dy;
	dy_ov_dt = 1.0/dt_ov_dy;
    
    // number of nodes of the primal and dual grid in the x-direction
	nx_p = params->n_space[0]+1;
	nx_d = params->n_space[0]+2;
    
    // number of nodes of the primal and dual grid in the y-direction
	ny_p = params->n_space[1]+1;
	ny_d = params->n_space[1]+2;
    
    
    // -----------------------------------------------------
	// Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    // West boundary
    double theta  = 0.0; //! \todo Introduce in parameters for Boundary cond., e.g., params->EMBoundary->theta_W
    double factor = 1.0 / (cos(theta) + dt_ov_dx);
    Alpha_SM_W    = 2.0                     * factor;
    Beta_SM_W     = (dt_ov_dx - cos(theta)) * factor;
    Gamma_SM_W    = 4.0 * cos(theta)        * factor;
    Delta_SM_W    = - (sin(theta)+dt_ov_dy) * factor;
    Epsilon_SM_W  = - (sin(theta)-dt_ov_dy) * factor;
    
    // East boundary
    theta        = M_PI;
    factor       = 1.0 / (cos(theta) - dt_ov_dx);
    Alpha_SM_E   = 2.0                     * factor;
    Beta_SM_E    = (dt_ov_dx + cos(theta)) * factor;
    Gamma_SM_E   = 4.0 * cos(theta)        * factor;
    Delta_SM_E   = (sin(theta)+dt_ov_dy)   * factor;
    Epsilon_SM_E = (sin(theta)-dt_ov_dy)   * factor;
    
    
    // ----------------------
	// Electromagnetic fields
	// ----------------------
    
	dimPrim.resize( params->nDim_field );
	dimDual.resize( params->nDim_field );
    
	// Dimension of the primal and dual grids
	for (size_t i=0 ; i<params->nDim_field ; i++) {
		// Standard scheme
		dimPrim[i] = params->n_space[i]+1;
		dimDual[i] = params->n_space[i]+2;
		// + Ghost domain
		dimPrim[i] += 2*params->oversize[i];
		dimDual[i] += 2*params->oversize[i];
	}
    
    // Allocation of the EM fields
	ostringstream name;
	name.str(""); name << "fex." << process_coord_x;
	Ex_  = new Field2D( dimPrim, 0, false, name.str() );
	name.str(""); name << "fey." << process_coord_x;
	Ey_  = new Field2D( dimPrim, 1, false, name.str() );
	name.str(""); name << "fez." << process_coord_x;
	Ez_  = new Field2D( dimPrim, 2, false, name.str() );
	name.str(""); name << "fbx." << process_coord_x;
	Bx_  = new Field2D( dimPrim, 0, true, name.str() );
	name.str(""); name << "fby." << process_coord_x;
	By_  = new Field2D( dimPrim, 1, true, name.str() );
	name.str(""); name << "fbz." << process_coord_x;
	Bz_  = new Field2D( dimPrim, 2, true, name.str() );
	Bx_m = new Field2D(dimPrim, 0, true);
	By_m = new Field2D(dimPrim, 1, true);
	Bz_m = new Field2D(dimPrim, 2, true);
	
	// Allocation of the total charge and currents
	name.str(""); name << "fjx." << process_coord_x;
	Jx_ = new Field2D( dimPrim, 0, false, name.str() );
	name.str(""); name << "fjy." << process_coord_x;
	Jy_ = new Field2D( dimPrim, 1, false, name.str() );
	name.str(""); name << "fjz." << process_coord_x;
	Jz_ = new Field2D( dimPrim, 2, false, name.str() );
	name.str(""); name << "rho." << process_coord_x;
	rho_ = new Field2D( dimPrim, name.str() );
	name.str(""); name << "rho_old." << process_coord_x;
	rho_o = new Field2D( dimPrim );
    
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    
    //!\ CHECK THIS: need to define for primal & dual grids (MG to JD)
	index_bc_min.resize( params->nDim_field, 0 );
	index_bc_max.resize( params->nDim_field, 0 );
	for (size_t i=0 ; i<params->nDim_field ; i++) {
		index_bc_min[i] = params->oversize[i];
		index_bc_max[i] = dimDual[i]-params->oversize[i]-1;
    }
    
} // END constructor Electromagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::~ElectroMagn2D()
{
    DEBUG( "to be implemented" );
    
}//END ElectroMagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Solve Poisson
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solvePoisson(SmileiMPI* smpi)
{
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
	Field2D* Ey2D = static_cast<Field2D*>(Ey_);
	Field2D* Ez2D = static_cast<Field2D*>(Ez_);
	Field2D* Bx2D = static_cast<Field2D*>(Bx_);
	Field2D* By2D = static_cast<Field2D*>(By_);
	Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    
    // AT THE MOMENT PUT ALL FIELDS TO 0 WHEN THEY ARE CREATED !!! (useful to test the Maxwell Solver)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D)(i,j) = 0.0;
        }
    }
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D)(i,j) = 0.0;
        }
    }
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D)(i,j) = 0.0;
        }
    }
    
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D)(i,j) = 0.0;
        }
    }
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(i,j) = 0.0;
        }
    }
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D)(i,j) = 0.0;
        }
    }
    
}//END solvePoisson


// ---------------------------------------------------------------------------------------------------------------------
// Maxwell solver using the FDTD scheme
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwell(double time_dual, SmileiMPI* smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    
    //! \todo All this is generic (does not depend on geometry) move to ElectroMagn.cpp ??? (MG to JD)
    saveMagneticFields();
    solveMaxwellAmpere();
    solveMaxwellFaraday();
    smpi2D->exchangeB( this );
    // ..
    // (*By1D)(0)       = (*By1D_west_neighbor)(dimDual-2*oversize)
    // (*By1D)(dimDual) = (*By1D_east_neighbor)(2*oversize)
    // ....
    applyEMBoundaryConditions(time_dual, smpi2D);
    centerMagneticFields();
    
}//END ElectroMagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::saveMagneticFields()
{
    // Static cast of the fields
	Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
	Field2D* By2D   = static_cast<Field2D*>(By_);
	Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
	Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
	Field2D* By2D_m = static_cast<Field2D*>(By_m);
	Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d)
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_m)(i,j)=(*Bx2D)(i,j);
        }
	}
    
    // Magnetic field By^(d,p)
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j)=(*By2D)(i,j);
        }
	}
	
    // Magnetic field Bz^(d,d)
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j)=(*Bz2D)(i,j);
        }
	}
    
}//END saveMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Ampere equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwellAmpere()
{
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
	Field2D* Ey2D = static_cast<Field2D*>(Ey_);
	Field2D* Ez2D = static_cast<Field2D*>(Ez_);
	Field2D* Bx2D = static_cast<Field2D*>(Bx_);
	Field2D* By2D = static_cast<Field2D*>(By_);
	Field2D* Bz2D = static_cast<Field2D*>(Bz_);
	Field2D* Jx2D = static_cast<Field2D*>(Jx_);
	Field2D* Jy2D = static_cast<Field2D*>(Jy_);
	Field2D* Jz2D = static_cast<Field2D*>(Jz_);
    
    // Electric field Ex^(d,p)
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D)(i,j) += -dt*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ey^(p,d)
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D)(i,j) += -dt*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ez^(p,p)
	for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D)(i,j) += -dt*(*Jz2D)(i,j)
            +               dt_ov_dx * ( (*By2D)(i+1,j) - (*By2D)(i,j) )
            -               dt_ov_dy * ( (*Bx2D)(i,j+1) - (*Bx2D)(i,j) );
        }
    }
    
}//END solveMaxwellAmpere



// ---------------------------------------------------------------------------------------------------------------------
// Solve the Maxwell-Faraday equation
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::solveMaxwellFaraday()
{
        
    // Static-cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
	Field2D* Ey2D = static_cast<Field2D*>(Ey_);
	Field2D* Ez2D = static_cast<Field2D*>(Ez_);
	Field2D* Bx2D = static_cast<Field2D*>(Bx_);
	Field2D* By2D = static_cast<Field2D*>(By_);
	Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    
	// Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p;  i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bx2D)(i,j) -= dt_ov_dy * ( (*Ez2D)(i,j) - (*Ez2D)(i,j-1) );
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(i,j) += dt_ov_dx * ( (*Ez2D)(i,j) - (*Ez2D)(i-1,j) );
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*Bz2D)(i,j) += dt_ov_dy * ( (*Ex2D)(i,j) - (*Ex2D)(i,j-1) )
            -               dt_ov_dx * ( (*Ey2D)(i,j) - (*Ey2D)(i-1,j) );
        }
    }
    
}//END solveMaxwellFaraday



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::applyEMBoundaryConditions(double time_dual, SmileiMPI* smpi)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    
    // Static cast of the fields
    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
	Field2D* Ey2D = static_cast<Field2D*>(Ey_);
	Field2D* Ez2D = static_cast<Field2D*>(Ez_);
	Field2D* Bx2D = static_cast<Field2D*>(Bx_);
	Field2D* By2D = static_cast<Field2D*>(By_);
	Field2D* Bz2D = static_cast<Field2D*>(Bz_);
    
    
    // ! \todo Transverse profile & incidence angle is not yet introduced (MG)
    // -----------------------------------------
    // Laser temporal profile
    // -----------------------------------------
	double byW=0.0, bzW=0.0, byE=0.0, bzE=0.0;
    
	for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
		// testing the time-profile
		// ------------------------
        
		if (laser_[ilaser]->laser_struct.time_profile == "constant") {
			if (laser_[ilaser]->laser_struct.angle == 0){
				// Incident field (left boundary)
				byW += laser_[ilaser]->a0_delta_y_ * sin(time_dual) ;//* laser_[ilaser]->time_profile(time_dual);
				bzW += laser_[ilaser]->a0_delta_z_ * cos(time_dual) ;//* laser_[ilaser]->time_profile(time_dual);
			} else if (laser_[ilaser]->laser_struct.angle == 180){
				// Incident field (right boundary)
				byE += laser_[ilaser]->a0_delta_y_ * sin(time_dual) ;//* laser_[ilaser]->time_profile(time_dual);
				bzE += laser_[ilaser]->a0_delta_z_ * cos(time_dual) ;//* laser_[ilaser]->time_profile(time_dual);
			} else {
				ERROR("Angle not yet implemented for laser " << ilaser);
			}
		} else {
			ERROR("Laser profile "<< ilaser <<" not allowed");
		}//ENDif time_profile
	}//ilaser
    
    
    // -----------------------------------------
    // Silver-Mueller boundary conditions (West)
    // -----------------------------------------
	if ( smpi2D->isWester() ) {

		// for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(index_bc_min[0],j) = Alpha_SM_W   * (*Ez2D)(index_bc_min[0],j)     // CHECK ALL i INDEXES !!!
            +                            Beta_SM_W    * (*By2D)(index_bc_min[0]+1,j)
            +                            Gamma_SM_W   * byW
            +                            Delta_SM_W   * (*Bx2D)(index_bc_max[0],j+1)
            +                            Epsilon_SM_W * (*Bx2D)(index_bc_max[0],j);
        }
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D)(index_bc_min[0],j) = -Alpha_SM_W * (*Ey2D)(index_bc_min[0],j)
            +                             Beta_SM_W  * (*Bz2D)(index_bc_min[0]+1,j)
            +                             Gamma_SM_W * bzW;
        }
        
		// Correction on unused extreme ghost cells: put the fields to 0
        // -------------------------------------------------------------
        // !\todo Take care that the fields are defined on different p/d grids ! (MG to JD)
        
        // for By^(d,p)
		for (unsigned int i=0 ; i<index_bc_min[0] ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*By2D)(i,j)=0.0;
            }
        }
        // for Bz^(d,d)
        for (unsigned int i=0 ; i<index_bc_min[0] ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bz2D)(i,j)=0.0;
            }
        }
        // for Ey^(p,d)
        for (unsigned int i=0 ; i<index_bc_min[0] ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Ey2D)(i,j)=0.0;
            }
        }
        // for Ez^(p,p)
        for (unsigned int i=0 ; i<index_bc_min[0] ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ez2D)(i,j)=0.0;
            }
		}
	}//if West
    
    
    // -----------------------------------------
    // Silver-Mueller boundary conditions (East)
    // -----------------------------------------
	if ( smpi2D->isEaster() ) {
        
		// for By^(d,p)
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D)(index_bc_max[0],j) = Alpha_SM_E   * (*Ez2D)(index_bc_max[0],j)   // CHECK, ALL i INDEXES (here, is it not index_bc_max[0] - 1 ???)
            +                            Beta_SM_E    * (*By2D)(index_bc_max[0]-1,j)
            +                            Gamma_SM_E   * byE
            +                            Delta_SM_E   * (*Bx2D)(index_bc_max[0]-1,j+1)
            +                            Epsilon_SM_E * (*Bx2D)(index_bc_max[0]-1,j);
        }
         
        // for Bz^(d,d)
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D)(index_bc_max[0],j) = -Alpha_SM_E * (*Ey2D)(index_bc_max[0],j)
            +                           Beta_SM_E    * (*Bz2D)(index_bc_max[0]-1,j)
            +                           Gamma_SM_E   * bzE;
        }
        
		// Correction on unused extreme ghost
        // ----------------------------------
        
        // for By^(d,p)
		for (unsigned int i=index_bc_max[0]+1 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*By2D)(i,j) = 0.0;
            }
        }
        // for Bz^(d,d)
		for (unsigned int i=index_bc_max[0]+1 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bz2D)(i,j) = 0.0;
            }
        }
        // for Ey^(p,d)
		for (unsigned int i=index_bc_max[0]+1 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Ey2D)(i,j)=0.0;
            }
		}
        // for Ez^(p,p)
		for (unsigned int i=index_bc_max[0]+1 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Ez2D)(i,j)=0.0;
            }
		}
	}//if East
    
}// END applyEMBoundaryConditions



// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::centerMagneticFields()
{
    // Static cast of the fields
	Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
	Field2D* By2D   = static_cast<Field2D*>(By_);
	Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
	Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
	Field2D* By2D_m = static_cast<Field2D*>(By_m);
	Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d)
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_m)(i,j) = ( (*Bx2D)(i,j) + (*Bx2D_m)(i,j) )*0.5;
        }
	}
    
    // Magnetic field By^(d,p)
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j) = ( (*By2D)(i,j) + (*By2D_m)(i,j) )*0.5;
        }
	}
	
    // Magnetic field Bz^(d,d)
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j) = ( (*Bz2D)(i,j) + (*Bz2D_m)(i,j) )*0.5;
        }
	}
    
}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::initRhoJ()
{
    
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
	Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
	Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
	Field2D* rho2D   = static_cast<Field2D*>(rho_);
	Field2D* rho2D_o = static_cast<Field2D*>(rho_o);
	
    // Save current charge density as old charge density
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*rho2D_o)(i,j) = (*rho2D)(i,j);
        }
	}
    
    // Charge density rho^(p,p) to 0
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*rho2D)(i,j) = 0.0;
        }
	}
    
    // Current Jx^(d,p) to 0
	for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jx2D)(i,j) = 0.0;
        }
	}
    
    // Current Jy^(p,d) to 0
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Jy2D)(i,j) = 0.0;
        }
	}
    
    // Current Jz^(p,p) to 0
	for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Jz2D)(i,j) = 0.0;
        }
	}
    
}//END initRhoJ
