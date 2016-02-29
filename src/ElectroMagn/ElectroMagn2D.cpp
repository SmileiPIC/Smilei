#include "ElectroMagn2D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "Params.h"
#include "Field2D.h"

#include "Patch.h"
#include <cstring>

#include "Profile.h"

#include "ElectroMagnBC.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::ElectroMagn2D(Params &params, vector<Species*>& vecSpecies, Patch* patch) : 
  ElectroMagn(params, vecSpecies, patch),
isWestern(patch->isWestern()),
isEastern(patch->isEastern()),
isNorthern(patch->isNorthern()),
isSouthern(patch->isSouthern())
{    
    
    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = cell_length[0];
    dt_ov_dx = timestep/dx;
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dy       = cell_length[1];
    dt_ov_dy = timestep/dy;
    dy_ov_dt = 1.0/dt_ov_dy;
    
    // ----------------------
    // Electromagnetic fields
    // ----------------------
    //! \todo Homogenize 1D/2D dimPrim/dimDual or nx_p/nx_d/ny_p/ny_d
    
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    
    // Dimension of the primal and dual grids
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1];
    
    // Allocation of the EM fields

    Ex_  = new Field2D(dimPrim, 0, false, "Ex");
    Ey_  = new Field2D(dimPrim, 1, false, "Ey");
    Ez_  = new Field2D(dimPrim, 2, false, "Ez");
    Bx_  = new Field2D(dimPrim, 0, true,  "Bx");
    By_  = new Field2D(dimPrim, 1, true,  "By");
    Bz_  = new Field2D(dimPrim, 2, true,  "Bz");
    Bx_m = new Field2D(dimPrim, 0, true,  "Bx_m");
    By_m = new Field2D(dimPrim, 1, true,  "By_m");
    Bz_m = new Field2D(dimPrim, 2, true,  "Bz_m");
    
    // Total charge currents and densities
    Jx_   = new Field2D(dimPrim, 0, false, "Jx");
    Jy_   = new Field2D(dimPrim, 1, false, "Jy");
    Jz_   = new Field2D(dimPrim, 2, false, "Jz");
    rho_  = new Field2D(dimPrim, "Rho" );
    
    // Allocation of the time-averaged EM fields
    Ex_avg  = new Field2D(dimPrim, 0, false, "Ex_avg");
    Ey_avg  = new Field2D(dimPrim, 1, false, "Ey_avg");
    Ez_avg  = new Field2D(dimPrim, 2, false, "Ez_avg");
    Bx_avg  = new Field2D(dimPrim, 0, true,  "Bx_avg");
    By_avg  = new Field2D(dimPrim, 1, true,  "By_avg");
    Bz_avg  = new Field2D(dimPrim, 2, true,  "Bz_avg");
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field2D(dimPrim, 0, false, ("Jx_"+vecSpecies[ispec]->species_type).c_str());
        Jy_s[ispec]  = new Field2D(dimPrim, 1, false, ("Jy_"+vecSpecies[ispec]->species_type).c_str());
        Jz_s[ispec]  = new Field2D(dimPrim, 2, false, ("Jz_"+vecSpecies[ispec]->species_type).c_str());
        rho_s[ispec] = new Field2D(dimPrim, ("Rho_"+vecSpecies[ispec]->species_type).c_str());
    }

    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (unsigned int i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    /*
     MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
            << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
     */
    
    
    // Define limits of non duplicated elements
    // (by construction 1 (prim) or 2 (dual) elements shared between per MPI process)
    // istart 
	for (unsigned int i=0 ; i<3 ; i++)
	    for (unsigned int isDual=0 ; isDual<2 ; isDual++)
		istart[i][isDual] = 0;
	for (unsigned int i=0 ; i<nDim_field ; i++) {
	    for (unsigned int isDual=0 ; isDual<2 ; isDual++) {
		istart[i][isDual] = oversize[i];
		if (patch->Pcoordinates[i]!=0) istart[i][isDual]+=1;
	    }
	}
    
	// bufsize = nelements
	for (unsigned int i=0 ; i<3 ; i++) 
	    for (unsigned int isDual=0 ; isDual<2 ; isDual++)
		bufsize[i][isDual] = 1;
    
	for (unsigned int i=0 ; i<nDim_field ; i++) {
	    for (int isDual=0 ; isDual<2 ; isDual++)
		bufsize[i][isDual] = n_space[i] + 1;
        
	    for (int isDual=0 ; isDual<2 ; isDual++) {
		bufsize[i][isDual] += isDual; 
		if ( params.number_of_patches[i]!=1 ) {                

		    if ( ( !isDual ) && (patch->Pcoordinates[i]!=0) )
			bufsize[i][isDual]--;
		    else if  (isDual) {
			bufsize[i][isDual]--;
			if ( (patch->Pcoordinates[i]!=0) && (patch->Pcoordinates[i]!=params.number_of_patches[i]-1) ) 
			    bufsize[i][isDual]--;
		    }
                
		} // if ( params.number_of_patches[i]!=1 )
	    } // for (int isDual=0 ; isDual
	} // for (unsigned int i=0 ; i<nDim_field 
    
    for (unsigned int i=0; i<antennas.size(); i++) {
        if      (antennas[i].fieldName == "Jx")
            antennas[i].field = new Field2D(dimPrim, 0, false, "Jx");
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = new Field2D(dimPrim, 1, false, "Jy");
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = new Field2D(dimPrim, 2, false, "Jz");
        
        if (antennas[i].field)
            applyExternalField(antennas[i].field, antennas[i].space_profile, patch);
    }
    
}//END constructor Electromagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn2D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn2D::~ElectroMagn2D()
{
}//END ElectroMagn2D



// ---------------------------------------------------------------------------------------------------------------------
// Begin of Solve Poisson methods
// ---------------------------------------------------------------------------------------------------------------------
// in VectorPatch::solvePoisson
//     - initPoisson
//     - compute_r
//     - compute_Ap
//     - compute_pAp
//     - update_pand_r
//     - update_p
//     - initE
//     - centeringE


void ElectroMagn2D::initPoisson(Patch *patch)
{
    Field2D* rho2D = static_cast<Field2D*>(rho_);

    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------

    index_min_p_.resize(2,0);
    index_max_p_.resize(2,0);
    
    index_min_p_[0] = oversize[0];
    index_min_p_[1] = oversize[1];
    index_max_p_[0] = nx_p - 2 - oversize[0];
    index_max_p_[1] = ny_p - 2 - oversize[1];
    if (patch->isWestern()) {
        index_min_p_[0] = 0;
    }
    if (patch->isEastern()) {
        index_max_p_[0] = nx_p-1;
    }

    phi_ = new Field2D(dimPrim);    // scalar potential
    r_   = new Field2D(dimPrim);    // residual vector
    p_   = new Field2D(dimPrim);    // direction vector
    Ap_  = new Field2D(dimPrim);    // A*p vector

    
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*phi_)(i,j)   = 0.0;
            (*r_)(i,j)     = -(*rho2D)(i,j);
            (*p_)(i,j)     = (*r_)(i,j);
        }//j
    }//i

} // initPoisson

double ElectroMagn2D::compute_r()
{
    double rnew_dot_rnew_local(0.);
    for (unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++) {
        for (unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++) {
            rnew_dot_rnew_local += (*r_)(i,j)*(*r_)(i,j);
        }
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn2D::compute_Ap(Patch* patch)
{
    double one_ov_dx_sq       = 1.0/(dx*dx);
    double one_ov_dy_sq       = 1.0/(dy*dy);
    double two_ov_dx2dy2      = 2.0*(1.0/(dx*dx)+1.0/(dy*dy));
    
    // vector product Ap = A*p
    for (unsigned int i=1; i<nx_p-1; i++) {
	for (unsigned int j=1; j<ny_p-1; j++) {
	    (*Ap_)(i,j) = one_ov_dx_sq*((*p_)(i-1,j)+(*p_)(i+1,j))
		+ one_ov_dy_sq*((*p_)(i,j-1)+(*p_)(i,j+1))
		- two_ov_dx2dy2*(*p_)(i,j);
	}//j
    }//i
        
        
    // Western BC
    if ( patch->isWestern() ) {
	for (unsigned int j=1; j<ny_p-1; j++) {
	    //Ap_(0,j)      = one_ov_dx_sq*(pWest[j]+p_(1,j))
	    (*Ap_)(0,j)      = one_ov_dx_sq*((*p_)(1,j))
                +              one_ov_dy_sq*((*p_)(0,j-1)+(*p_)(0,j+1))
                -              two_ov_dx2dy2*(*p_)(0,j);
	}
	// at corners
	//Ap_(0,0)           = one_ov_dx_sq*(pWest[0]+p_(1,0))               // West/South
        //    +                   one_ov_dy_sq*(pSouth[0]+p_(0,1))
	(*Ap_)(0,0)           = one_ov_dx_sq*((*p_)(1,0))               // West/South
            +                   one_ov_dy_sq*((*p_)(0,1))
            -                   two_ov_dx2dy2*(*p_)(0,0);
	//Ap_(0,ny_p-1)      = one_ov_dx_sq*(pWest[ny_p-1]+p_(1,ny_p-1))     // West/North
        //    +                   one_ov_dy_sq*(p_(0,ny_p-2)+pNorth[0])
	(*Ap_)(0,ny_p-1)      = one_ov_dx_sq*((*p_)(1,ny_p-1))     // West/North
            +                   one_ov_dy_sq*((*p_)(0,ny_p-2))
            -                   two_ov_dx2dy2*(*p_)(0,ny_p-1);
    }
        
    // Eastern BC
    if ( patch->isEastern() ) {
            
	for (unsigned int j=1; j<ny_p-1; j++) {
	    //Ap_(nx_p-1,j) = one_ov_dx_sq*(p_(nx_p-2,j)+pEast[j])
	    (*Ap_)(nx_p-1,j) = one_ov_dx_sq*((*p_)(nx_p-2,j))
                +              one_ov_dy_sq*((*p_)(nx_p-1,j-1)+(*p_)(nx_p-1,j+1))
                -              two_ov_dx2dy2*(*p_)(nx_p-1,j);
	}
	// at corners
	//Ap_(nx_p-1,0)      = one_ov_dx_sq*(p_(nx_p-2,0)+pEast[0])                 // East/South
        //    +                   one_ov_dy_sq*(pSouth[nx_p-1]+p_(nx_p-1,1))
	(*Ap_)(nx_p-1,0)      = one_ov_dx_sq*((*p_)(nx_p-2,0))                 // East/South
            +                   one_ov_dy_sq*((*p_)(nx_p-1,1))
            -                   two_ov_dx2dy2*(*p_)(nx_p-1,0);
	//Ap_(nx_p-1,ny_p-1) = one_ov_dx_sq*(p_(nx_p-2,ny_p-1)+pEast[ny_p-1])       // East/North
        //    +                   one_ov_dy_sq*(p_(nx_p-1,ny_p-2)+pNorth[nx_p-1])
	(*Ap_)(nx_p-1,ny_p-1) = one_ov_dx_sq*((*p_)(nx_p-2,ny_p-1))       // East/North
            +                   one_ov_dy_sq*((*p_)(nx_p-1,ny_p-2))
            -                   two_ov_dx2dy2*(*p_)(nx_p-1,ny_p-1);
    }
        
} // compute_pAp

double ElectroMagn2D::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    for (unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++) {
	for (unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++) {
	    p_dot_Ap_local += (*p_)(i,j)*(*Ap_)(i,j);
	}
    }
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagn2D::update_pand_r(double r_dot_r, double p_dot_Ap)
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for(unsigned int i=0; i<nx_p; i++) {
	for(unsigned int j=0; j<ny_p; j++) {
	    (*phi_)(i,j) += alpha_k * (*p_)(i,j);
	    (*r_)(i,j)   -= alpha_k * (*Ap_)(i,j);
	}
    }

} // update_pand_r

void ElectroMagn2D::update_p(double rnew_dot_rnew, double r_dot_r)
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for (unsigned int i=0; i<nx_p; i++) {
	for(unsigned int j=0; j<ny_p; j++) {
	    (*p_)(i,j) = (*r_)(i,j) + beta_k * (*p_)(i,j);
	}
    }
} // update_p

void ElectroMagn2D::initE(Patch *patch)
{
    Field2D* Ex2D  = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D  = static_cast<Field2D*>(Ey_);
    Field2D* rho2D = static_cast<Field2D*>(rho_);

    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    
    // Ex
    DEBUG("Computing Ex from scalar potential");
    for (unsigned int i=1; i<nx_d-1; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(i,j) = ((*phi_)(i-1,j)-(*phi_)(i,j))/dx;
        }
    }
    // Ey
    DEBUG("Computing Ey from scalar potential");
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=1; j<ny_d-1; j++) {
            (*Ey2D)(i,j) = ((*phi_)(i,j-1)-(*phi_)(i,j))/dy;
        }
    }

    // Apply BC on Ex and Ey
    // ---------------------
    // Ex / West
    if (patch->isWestern()) {
        DEBUG("Computing Western BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(0,j) = (*Ex2D)(1,j) + ((*Ey2D)(0,j+1)-(*Ey2D)(0,j))*dx/dy  - dx*(*rho2D)(0,j);
        }
    }
    // Ex / East
    if (patch->isEastern()) {
        DEBUG("Computing Eastern BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(nx_d-1,j) = (*Ex2D)(nx_d-2,j) - ((*Ey2D)(nx_p-1,j+1)-(*Ey2D)(nx_p-1,j))*dx/dy + dx*(*rho2D)(nx_p-1,j);
        }
    }

    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;

} // initE


void ElectroMagn2D::centeringE( std::vector<double> E_Add )
{
    Field2D* Ex2D  = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D  = static_cast<Field2D*>(Ey_);

    // Centering electrostatic fields
    for (unsigned int i=0; i<nx_d; i++) {
	for (unsigned int j=0; j<ny_p; j++) {
	    (*Ex2D)(i,j) += E_Add[0];
	}
    }
    for (unsigned int i=0; i<nx_p; i++) {
	for (unsigned int j=0; j<ny_d; j++) {
	    (*Ey2D)(i,j) += E_Add[1];
	}
    }
} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods 
// ---------------------------------------------------------------------------------------------------------------------


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
        memcpy(&((*Bx2D_m)(i,0)), &((*Bx2D)(i,0)),ny_d*sizeof(double) );
        //for (unsigned int j=0 ; j<ny_d ; j++) {
        //    (*Bx2D_m)(i,j)=(*Bx2D)(i,j);
        //}
    
    // Magnetic field By^(d,p)
        memcpy(&((*By2D_m)(i,0)), &((*By2D)(i,0)),ny_p*sizeof(double) );
        //for (unsigned int j=0 ; j<ny_p ; j++) {
        //    (*By2D_m)(i,j)=(*By2D)(i,j);
        //}
    
    // Magnetic field Bz^(d,d)
        memcpy(&((*Bz2D_m)(i,0)), &((*Bz2D)(i,0)),ny_d*sizeof(double) );
        //for (unsigned int j=0 ; j<ny_d ; j++) {
        //    (*Bz2D_m)(i,j)=(*Bz2D)(i,j);
        //}
    }// end for i
        memcpy(&((*By2D_m)(nx_p,0)), &((*By2D)(nx_p,0)),ny_p*sizeof(double) );
        //for (unsigned int j=0 ; j<ny_p ; j++) {
        //    (*By2D_m)(nx_p,j)=(*By2D)(nx_p,j);
        //}
        memcpy(&((*Bz2D_m)(nx_p,0)), &((*Bz2D)(nx_p,0)),ny_d*sizeof(double) );
        //for (unsigned int j=0 ; j<ny_d ; j++) {
        //    (*Bz2D_m)(nx_p,j)=(*Bz2D)(nx_p,j);
        //}
    
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
            (*Ex2D)(i,j) += -timestep*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ey^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D)(i,j) += -timestep*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) );
        }
    }
    
    // Electric field Ez^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D)(i,j) += -timestep*(*Jz2D)(i,j)
            +               dt_ov_dx * ( (*By2D)(i+1,j) - (*By2D)(i,j) )
            -               dt_ov_dy * ( (*Bx2D)(i,j+1) - (*Bx2D)(i,j) );
        }
    }
#ifdef _PATCH_DEBUG
    cout << "\tEx = "  << Ex_->norm() << endl;
    cout << "\tEy = "  << Ey_->norm() << endl;
#endif

}//END solveMaxwellAmpere


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
//    }
    
    // Magnetic field By^(d,p)
//    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j) = ( (*By2D)(i,j) + (*By2D_m)(i,j) )*0.5;
        }
//    }
    
    // Magnetic field Bz^(d,d)
//    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j) = ( (*Bz2D)(i,j) + (*Bz2D_m)(i,j) )*0.5;
        } // end for j
      } // end for i
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(nx_p,j) = ( (*By2D)(nx_p,j) + (*By2D_m)(nx_p,j) )*0.5;
        }
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(nx_p,j) = ( (*Bz2D)(nx_p,j) + (*Bz2D_m)(nx_p,j) )*0.5;
        } // end for j
/*    cout << "\tBx_m = "  << Bx_m->norm() << endl;
      cout << "\tBy_m = "  << By_m->norm() << endl;*/
#ifdef _PATCH_DEBUG
    cout << "\tBz_m = "  << Bz_m->norm() << endl;
#endif

    
}//END centerMagneticFields



// ---------------------------------------------------------------------------------------------------------------------
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg)
{
    // Static cast of the fields
    Field2D* Ex2D     = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D     = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D     = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D_m   = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m   = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m   = static_cast<Field2D*>(Bz_m);
    Field2D* Ex2D_avg = static_cast<Field2D*>(Ex_avg);
    Field2D* Ey2D_avg = static_cast<Field2D*>(Ey_avg);
    Field2D* Ez2D_avg = static_cast<Field2D*>(Ez_avg);
    Field2D* Bx2D_avg = static_cast<Field2D*>(Bx_avg);
    Field2D* By2D_avg = static_cast<Field2D*>(By_avg);
    Field2D* Bz2D_avg = static_cast<Field2D*>(Bz_avg);
    
    // reset the averaged fields for (time_step-1)%ntime_step_avg == 0
    if ( (time_step-1)%ntime_step_avg==0 ){
        Ex2D_avg->put_to(0.0);
        Ey2D_avg->put_to(0.0);
        Ez2D_avg->put_to(0.0);
        Bx2D_avg->put_to(0.0);
        By2D_avg->put_to(0.0);
        Bz2D_avg->put_to(0.0);
    }
    
    // increment the time-averaged fields
    
    // Electric field Ex^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ex2D_avg)(i,j) += (*Ex2D)(i,j);
        }
    }
    
    // Electric field Ey^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Ey2D_avg)(i,j) += (*Ey2D)(i,j);
        }
    }
    
    // Electric field Ez^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*Ez2D_avg)(i,j) += (*Ez2D)(i,j);
        }
    }
    
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_avg)(i,j) += (*Bx2D_m)(i,j);
        }
    }
    
    // Magnetic field By^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_avg)(i,j) += (*By2D_m)(i,j);
        }
    }
    
    // Magnetic field Bz^(d,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_avg)(i,j) += (*Bz2D_m)(i,j);
        }
    }
    
    
}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Reinitialize the total charge densities and currents
// - save current density as old density (charge conserving scheme)
// - put the new density and currents to 0
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::restartRhoJ()
{
    // --------------------------
    // Total currents and density
    // --------------------------
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);
        
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
    
}//END restartRhoJ
    
    
void ElectroMagn2D::restartRhoJs()
{
    unsigned int size_proj_buffer = 2*oversize[0]+clrw + 1 ;

    for (unsigned int ispec=0 ; ispec < n_species ; ispec++) {
        // -----------------------------------
        // Species currents and charge density
        // -----------------------------------
        Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
        Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
        Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
        Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);
   
        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*rho2D_s)(i,j) = 0.0;
            }
        }
        // Current Jx^(d,p) to 0
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jx2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jy^(p,d) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Jy2D_s)(i,j) = 0.0;
            }
        }
        
        // Current Jz^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Jz2D_s)(i,j) = 0.0;
            }
        }
    }//End loop on species.
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);
    
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

}//END restartRhoJs
    



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::computeTotalRhoJ()
{
    // static cast of the total currents and densities
    Field2D* Jx2D    = static_cast<Field2D*>(Jx_);
    Field2D* Jy2D    = static_cast<Field2D*>(Jy_);
    Field2D* Jz2D    = static_cast<Field2D*>(Jz_);
    Field2D* rho2D   = static_cast<Field2D*>(rho_);
    
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
        Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
        Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
        Field2D* rho2D_s = static_cast<Field2D*>(rho_s[ispec]);
        
        // Charge density rho^(p,p) to 0
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                if(ispec==0){
                    (*rho2D)(i,j) = 0.;
                    (*Jx2D)(i,j)  = 0.;
                    (*Jy2D)(i,j)  = 0.;
                    (*Jz2D)(i,j)  = 0.;
                }
                (*rho2D)(i,j) += (*rho2D_s)(i,j);
                (*Jx2D)(i,j) += (*Jx2D_s)(i,j);
                (*Jy2D)(i,j) += (*Jy2D_s)(i,j);
                (*Jz2D)(i,j) += (*Jz2D_s)(i,j);
            }
            if(ispec==0) (*Jy2D)(i,ny_p) = 0. ;
            (*Jy2D)(i,ny_p) += (*Jy2D_s)(i,ny_p);
        }
        
        {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                if(ispec==0) (*Jx2D)(nx_p,j) = 0. ;
                (*Jx2D)(nx_p,j) += (*Jx2D_s)(nx_p,j);
            }
        }
        
    }//END loop on species ispec
//END computeTotalRhoJ
}


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn2D::computePoynting() {

    Field2D* Ex2D     = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D     = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D     = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D_m   = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m   = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m   = static_cast<Field2D*>(Bz_m);

    if (isWestern) {
        unsigned int iEy=istart[0][Ey2D->isDual(0)];
        unsigned int iBz=istart[0][Bz2D_m->isDual(0)];
        unsigned int iEz=istart[0][Ez2D->isDual(0)];
        unsigned int iBy=istart[0][By2D_m->isDual(0)];
        
        unsigned int jEy=istart[1][Ey2D->isDual(1)];
        unsigned int jBz=istart[1][Bz2D_m->isDual(1)];
        unsigned int jEz=istart[1][Ez2D->isDual(1)];
        unsigned int jBy=istart[1][By2D_m->isDual(1)];
        
        for (unsigned int j=0; j<=bufsize[1][Ez2D->isDual(1)]; j++) {
            
            double Ey__ = 0.5*((*Ey2D)(iEy,jEy+j) + (*Ey2D)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz2D_m)(iBz,jBz+j)+(*Bz2D_m)(iBz+1,jBz+j)+(*Bz2D_m)(iBz,jBz+j+1)+(*Bz2D_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez2D)(iEz,jEz+j);
            double By__ = 0.5*((*By2D_m)(iBy,jBy+j) + (*By2D_m)(iBy+1, jBy+j));
            
            poynting_inst[0][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[0][0]+= poynting_inst[0][0];
        }
    }//if Western
    
    
    if (isEastern) {
        unsigned int iEy=istart[0][Ey2D->isDual(0)]  + bufsize[0][Ey2D->isDual(0)] -1;
        unsigned int iBz=istart[0][Bz2D_m->isDual(0)] + bufsize[0][Bz2D_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez2D->isDual(0)]  + bufsize[0][Ez2D->isDual(0)] -1;
        unsigned int iBy=istart[0][By2D_m->isDual(0)] + bufsize[0][By2D_m->isDual(0)]-1;
        
        unsigned int jEy=istart[1][Ey2D->isDual(1)];
        unsigned int jBz=istart[1][Bz2D_m->isDual(1)];
        unsigned int jEz=istart[1][Ez2D->isDual(1)];
        unsigned int jBy=istart[1][By2D_m->isDual(1)];
        
        for (unsigned int j=0; j<=bufsize[1][Ez2D->isDual(1)]; j++) {
            
            double Ey__ = 0.5*((*Ey2D)(iEy,jEy+j) + (*Ey2D)(iEy, jEy+j+1));
            double Bz__ = 0.25*((*Bz2D_m)(iBz,jBz+j)+(*Bz2D_m)(iBz+1,jBz+j)+(*Bz2D_m)(iBz,jBz+j+1)+(*Bz2D_m)(iBz+1,jBz+j+1));
            double Ez__ = (*Ez2D)(iEz,jEz+j);
            double By__ = 0.5*((*By2D_m)(iBy,jBy+j) + (*By2D_m)(iBy+1, jBy+j));
            
            poynting_inst[1][0] = dy*timestep*(Ey__*Bz__ - Ez__*By__);
            poynting[1][0] -= poynting_inst[1][0];
        }
    }//if Easter
    
    if (isSouthern) {
        
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)]; 
        
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];
        
        for (unsigned int i=0; i<=bufsize[0][Ez2D->isDual(0)]; i++) {
            double Ez__ = (*Ez2D)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
            
            poynting_inst[0][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[0][1] += poynting_inst[0][1];
        }
    }// if South
    
    if (isNorthern) {
        unsigned int iEz=istart[0][Ez2D->isDual(0)];
        unsigned int iBx=istart[0][Bx2D_m->isDual(0)];
        unsigned int iEx=istart[0][Ex2D->isDual(0)];
        unsigned int iBz=istart[0][Bz2D_m->isDual(0)];
        
        unsigned int jEz=istart[1][Ez2D->isDual(1)]  + bufsize[1][Ez2D->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx2D_m->isDual(1)] + bufsize[1][Bx2D_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex2D->isDual(1)]  + bufsize[1][Ex2D->isDual(1)] -1;
        unsigned int jBz=istart[1][Bz2D_m->isDual(1)] + bufsize[1][Bz2D_m->isDual(1)]-1;
        
        for (unsigned int i=0; i<=bufsize[0][Ez_->isDual(0)]; i++) {
            double Ez__ = (*Ez2D)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
            
            poynting_inst[1][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[1][1] -= poynting_inst[1][1];
        }
    }//if North

}

void ElectroMagn2D::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    
    Field2D* field2D=static_cast<Field2D*>(my_field);
    
    vector<double> pos(2,0);
    pos[0] = (double)(patch->getCellStartingGlobalIndex(0)+(field2D->isDual(0)?-0.5:0));
    pos[1] = (double)(patch->getCellStartingGlobalIndex(1)+(field2D->isDual(1)?-0.5:0));
    int N0 = (int)field2D->dims()[0];
    int N1 = (int)field2D->dims()[1];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (int i=0 ; i<N0 ; i++) {
        pos[0] += dx;
        for (int j=0 ; j<N1 ; j++) {
            pos[1] += dy;
            (*field2D)(i,j) += profile->valueAt(pos);
        }
    }
    
    if (emBoundCond[0]!=0) emBoundCond[0]->save_fields_BC2D_Long(my_field);
    if (emBoundCond[1]!=0) emBoundCond[1]->save_fields_BC2D_Long(my_field);
    if (emBoundCond[2]!=0) emBoundCond[2]->save_fields_BC2D_Trans(my_field);
    if (emBoundCond[3]!=0) emBoundCond[3]->save_fields_BC2D_Trans(my_field);
    
}
