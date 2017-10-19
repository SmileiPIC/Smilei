#include "ElectroMagn3D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "Params.h"
#include "Field3D.h"

#include "Patch.h"
#include <cstring>

#include "Profile.h"

#include "ElectroMagnBC.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::ElectroMagn3D(Params &params, DomainDecomposition* domain_decomposition, vector<Species*>& vecSpecies, Patch* patch) : 
  ElectroMagn(params, domain_decomposition, vecSpecies, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax()),
isZmax(patch->isZmax()),
isZmin(patch->isZmin())
{    
    
    initElectroMagn3DQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field3D(("Jx_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        Jy_s[ispec]  = new Field3D(("Jy_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        Jz_s[ispec]  = new Field3D(("Jz_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        rho_s[ispec] = new Field3D(("Rho_"+vecSpecies[ispec]->species_type).c_str(), dimPrim);
    }
    
}//END constructor Electromagn3D


ElectroMagn3D::ElectroMagn3D( ElectroMagn3D* emFields, Params &params, Patch* patch ) : 
    ElectroMagn(emFields, params, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax()),
isZmax(patch->isZmax()),
isZmin(patch->isZmin())
{
    
    initElectroMagn3DQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        if ( emFields->Jx_s[ispec] != NULL ) {
            if ( emFields->Jx_s[ispec]->data_ != NULL )
                Jx_s[ispec]  = new Field3D(dimPrim, 0, false, emFields->Jx_s[ispec]->name);
            else
                Jx_s[ispec]  = new Field3D(emFields->Jx_s[ispec]->name, dimPrim);
        }
        if ( emFields->Jy_s[ispec] != NULL ) {
            if ( emFields->Jy_s[ispec]->data_ != NULL )
                Jy_s[ispec]  = new Field3D(dimPrim, 1, false, emFields->Jy_s[ispec]->name);
            else
                Jy_s[ispec]  = new Field3D(emFields->Jy_s[ispec]->name, dimPrim);
        }
        if ( emFields->Jz_s[ispec] != NULL ) {
            if ( emFields->Jz_s[ispec]->data_ != NULL )
                Jz_s[ispec]  = new Field3D(dimPrim, 2, false, emFields->Jz_s[ispec]->name);
            else
                Jz_s[ispec]  = new Field3D(emFields->Jz_s[ispec]->name, dimPrim);
        }
        if ( emFields->rho_s[ispec] != NULL ) {
            if ( emFields->rho_s[ispec]->data_ != NULL )
                rho_s[ispec] = new Field3D(dimPrim, emFields->rho_s[ispec]->name );
            else
                rho_s[ispec]  = new Field3D(emFields->rho_s[ispec]->name, dimPrim);
        }
    }


}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn3D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::initElectroMagn3DQuantities(Params &params, Patch* patch)
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
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the z-direction)
    dz       = cell_length[2];
    dt_ov_dz = timestep/dz;
    dz_ov_dt = 1.0/dt_ov_dz;
    
    // ----------------------
    // Electromagnetic fields
    // ----------------------
    
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
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = n_space[2]+1+2*oversize[2];
    nz_d = n_space[2]+2+2*oversize[2];
    
    // Allocation of the EM fields
    
    Ex_  = new Field3D(dimPrim, 0, false, "Ex");
    Ey_  = new Field3D(dimPrim, 1, false, "Ey");
    Ez_  = new Field3D(dimPrim, 2, false, "Ez");
    Bx_  = new Field3D(dimPrim, 0, true,  "Bx");
    By_  = new Field3D(dimPrim, 1, true,  "By");
    Bz_  = new Field3D(dimPrim, 2, true,  "Bz");
    Bx_m = new Field3D(dimPrim, 0, true,  "Bx_m");
    By_m = new Field3D(dimPrim, 1, true,  "By_m");
    Bz_m = new Field3D(dimPrim, 2, true,  "Bz_m");
    
    // Total charge currents and densities
    Jx_   = new Field3D(dimPrim, 0, false, "Jx");
    Jy_   = new Field3D(dimPrim, 1, false, "Jy");
    Jz_   = new Field3D(dimPrim, 2, false, "Jz");
    rho_  = new Field3D(dimPrim, "Rho" );
    if(params.is_pxr == true) {
        rhoold_ = new Field3D(dimPrim,"Rho");
        Ex_pxr  = new Field3D(dimDual);
        Ey_pxr  = new Field3D(dimDual);
        Ez_pxr  = new Field3D(dimDual);
        Bx_pxr  = new Field3D(dimDual);
        By_pxr  = new Field3D(dimDual);
        Bz_pxr  = new Field3D(dimDual);
        Jx_pxr  = new Field3D(dimDual);
        Jy_pxr  = new Field3D(dimDual);
        Jz_pxr  = new Field3D(dimDual);
        rho_pxr = new Field3D(dimDual);
        rhoold_pxr  = new Field3D(dimDual);
	std::cout<<"dimdual"<<dimDual[0]<<" "<<dimDual[1]<<" "<<dimDual[2]<<std::endl;

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
}

// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::~ElectroMagn3D()
{
}//END ElectroMagn3D



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


void ElectroMagn3D::initPoisson(Patch *patch)
{
    Field3D* rho3D = static_cast<Field3D*>(rho_);

    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------

    index_min_p_.resize(3,0);
    index_max_p_.resize(3,0);
    
    index_min_p_[0] = oversize[0];
    index_min_p_[1] = oversize[1];
    index_min_p_[2] = oversize[2];
    index_max_p_[0] = nx_p - 2 - oversize[0];
    index_max_p_[1] = ny_p - 2 - oversize[1];
    index_max_p_[2] = nz_p - 2 - oversize[2];
    if (patch->isXmin()) {
        index_min_p_[0] = 0;
    }
    if (patch->isXmax()) {
        index_max_p_[0] = nx_p-1;
    }

    phi_ = new Field3D(dimPrim);    // scalar potential
    r_   = new Field3D(dimPrim);    // residual vector
    p_   = new Field3D(dimPrim);    // direction vector
    Ap_  = new Field3D(dimPrim);    // A*p vector

    
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*phi_)(i,j,k)   = 0.0;
                (*r_)(i,j,k)     = -(*rho3D)(i,j,k);
                (*p_)(i,j,k)     = (*r_)(i,j,k);
            }
        }//j
    }//i

} // initPoisson

double ElectroMagn3D::compute_r()
{
    double rnew_dot_rnew_local(0.);
    for (unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++) {
        for (unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++) {
            for (unsigned int k=index_min_p_[2]; k<=index_max_p_[2]; k++) {
                rnew_dot_rnew_local += (*r_)(i,j,k)*(*r_)(i,j,k);
            }
        }
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn3D::compute_Ap(Patch* patch)
{
    double one_ov_dx_sq       = 1.0/(dx*dx);
    double one_ov_dy_sq       = 1.0/(dy*dy);
    double one_ov_dz_sq       = 1.0/(dz*dz);
    double three_ov_dx2dy2dz2 = 3.0*(1.0/(dx*dx)+1.0/(dy*dy)+1.0/(dz*dz));
    
    // vector product Ap = A*p
    for (unsigned int i=1; i<nx_p-1; i++) {
        for (unsigned int j=1; j<ny_p-1; j++) {
            for (unsigned int k=1; k<nz_p-1; k++) {
                (*Ap_)(i,j,k) = one_ov_dx_sq*((*p_)(i-1,j,k)+(*p_)(i+1,j,k))
                    + one_ov_dy_sq*((*p_)(i,j-1,k)+(*p_)(i,j+1,k))
                    + one_ov_dz_sq*((*p_)(i,j,k-1)+(*p_)(i,j,k+1))
                    - three_ov_dx2dy2dz2*(*p_)(i,j,k);
            }//k
        }//j
    }//i
        
        
    // Xmin BC
    if ( patch->isXmin() ) {
        for (unsigned int j=1; j<ny_p-1; j++) {
            for (unsigned int k=1; k<nz_p-1; k++) {
                (*Ap_)(0,j,k)      = one_ov_dx_sq*((*p_)(1,j,k))
                    +              one_ov_dy_sq*((*p_)(0,j-1,k)+(*p_)(0,j+1,k))
                    +              one_ov_dz_sq*((*p_)(0,j,k-1)+(*p_)(0,j,k+1))
                    -              three_ov_dx2dy2dz2*(*p_)(0,j,k);
            }
        }
        // at corners
        (*Ap_)(0,0,0)           = one_ov_dx_sq*((*p_)(1,0,0))         // Xmin/Ymin/Zmin
            +                   one_ov_dy_sq*((*p_)(0,1,0))
            +                   one_ov_dz_sq*((*p_)(0,0,1))
            -                   three_ov_dx2dy2dz2*(*p_)(0,0,0);
        (*Ap_)(0,ny_p-1,0)      = one_ov_dx_sq*((*p_)(1,ny_p-1,0))     // Xmin/Ymax/Zmin
            +                   one_ov_dy_sq*((*p_)(0,ny_p-2,0))
            +                   one_ov_dz_sq*((*p_)(0,ny_p-1,1))
            -                   three_ov_dx2dy2dz2*(*p_)(0,ny_p-1,0);
        (*Ap_)(0,0,nz_p-1)      = one_ov_dx_sq*((*p_)(1,0,nz_p-1))          // Xmin/Ymin/Zmin
            +                   one_ov_dy_sq*((*p_)(0,1,nz_p-1))
            +                   one_ov_dz_sq*((*p_)(0,0,nz_p-2))
            -                   three_ov_dx2dy2dz2*(*p_)(0,0,nz_p-1);
        (*Ap_)(0,ny_p-1,nz_p-1) = one_ov_dx_sq*((*p_)(1,ny_p-1,nz_p-1))     // Xmin/Ymax/Zmin
            +                   one_ov_dy_sq*((*p_)(0,ny_p-2,nz_p-1))
            +                   one_ov_dz_sq*((*p_)(0,ny_p-1,nz_p-2))
            -                   three_ov_dx2dy2dz2*(*p_)(0,ny_p-1,nz_p-1);
    }
        
    // Xmax BC
    if ( patch->isXmax() ) {
            
        for (unsigned int j=1; j<ny_p-1; j++) {
            for (unsigned int k=1; k<nz_p-1; k++) {
                (*Ap_)(nx_p-1,j,k) = one_ov_dx_sq*((*p_)(nx_p-2,j,k))
                    +              one_ov_dy_sq*((*p_)(nx_p-1,j-1,k)+(*p_)(nx_p-1,j+1,k))
                    +              one_ov_dz_sq*((*p_)(nx_p-1,j,k-1)+(*p_)(nx_p-1,j,k+1))
                    -              three_ov_dx2dy2dz2*(*p_)(nx_p-1,j,k);
            }
        }
        // at corners
        (*Ap_)(nx_p-1,0,0)      = one_ov_dx_sq*((*p_)(nx_p-2,0,0))            // Xmax/Ymin/Zmin
            +                   one_ov_dy_sq*((*p_)(nx_p-1,1,0))
            +                   one_ov_dz_sq*((*p_)(nx_p-1,0,1))
            -                   three_ov_dx2dy2dz2*(*p_)(nx_p-1,0,0);
        (*Ap_)(nx_p-1,ny_p-1,0) = one_ov_dx_sq*((*p_)(nx_p-2,ny_p-1,0))       // Xmax/Ymax/Zmin
            +                   one_ov_dy_sq*((*p_)(nx_p-1,ny_p-2,0))
            +                   one_ov_dz_sq*((*p_)(nx_p-1,ny_p-1,1))
            -                   three_ov_dx2dy2dz2*(*p_)(nx_p-1,ny_p-1,0);
        (*Ap_)(nx_p-1,0,nz_p-1)      = one_ov_dx_sq*((*p_)(nx_p-2,0,0))             // Xmax/Ymin/Zmax
            +                   one_ov_dy_sq*((*p_)(nx_p-1,1,nz_p-1))
            +                   one_ov_dz_sq*((*p_)(nx_p-1,0,nz_p-2))
            -                   three_ov_dx2dy2dz2*(*p_)(nx_p-1,0,nz_p-1);
        (*Ap_)(nx_p-1,ny_p-1,nz_p-1) = one_ov_dx_sq*((*p_)(nx_p-2,ny_p-1,nz_p-1))       // Xmax/Ymax/Zmax
            +                   one_ov_dy_sq*((*p_)(nx_p-1,ny_p-2,nz_p-1))
            +                   one_ov_dz_sq*((*p_)(nx_p-1,ny_p-1,nz_p-2))
            -                   three_ov_dx2dy2dz2*(*p_)(nx_p-1,ny_p-1,nz_p-1);
    }
        
} // compute_pAp

double ElectroMagn3D::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    for (unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++) {
        for (unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++) {
            for (unsigned int k=index_min_p_[2]; k<=index_max_p_[2]; k++) {
                p_dot_Ap_local += (*p_)(i,j,k)*(*Ap_)(i,j,k);
            }
        }
    }
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagn3D::update_pand_r(double r_dot_r, double p_dot_Ap)
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for(unsigned int i=0; i<nx_p; i++) {
        for(unsigned int j=0; j<ny_p; j++) {
            for(unsigned int k=0; k<nz_p; k++) {
                (*phi_)(i,j,k) += alpha_k * (*p_)(i,j,k);
                (*r_)(i,j,k)   -= alpha_k * (*Ap_)(i,j,k);
            }
        }
    }

} // update_pand_r

void ElectroMagn3D::update_p(double rnew_dot_rnew, double r_dot_r)
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for (unsigned int i=0; i<nx_p; i++) {
        for(unsigned int j=0; j<ny_p; j++) {
            for(unsigned int k=0; k<nz_p; k++) {
                (*p_)(i,j,k) = (*r_)(i,j,k) + beta_k * (*p_)(i,j,k);
            }
        }
    }
} // update_p

void ElectroMagn3D::initE(Patch *patch)
{
    Field3D* Ex3D  = static_cast<Field3D*>(Ex_);
    Field3D* Ey3D  = static_cast<Field3D*>(Ey_);
    Field3D* Ez3D  = static_cast<Field3D*>(Ez_);
    Field3D* rho3D = static_cast<Field3D*>(rho_);

    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    // Ex
    DEBUG("Computing Ex from scalar potential");
    for (unsigned int i=1; i<nx_d-1; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ex3D)(i,j,k) = ((*phi_)(i-1,j,k)-(*phi_)(i,j,k))/dx;
            }
        }
    }
    // Ey
    DEBUG("Computing Ey from scalar potential");
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=1; j<ny_d-1; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ey3D)(i,j,k) = ((*phi_)(i,j-1,k)-(*phi_)(i,j,k))/dy;
            }
        }
    }
    DEBUG("Computing Ez from scalar potential");
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=1; k<nz_d-1; k++) {
                (*Ez3D)(i,j,k) = ((*phi_)(i,j,k-1)-(*phi_)(i,j,k))/dz;
            }
        }
    }

    // Apply BC on Ex and Ey
    // ---------------------
    // Ex / Xmin
    if (patch->isXmin()) {
        DEBUG("Computing Xmin BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ex3D)(0,j,k) = (*Ex3D)(1,j,k) - dx*(*rho3D)(0,j,k)
                    + ((*Ey3D)(0,j+1,k)-(*Ey3D)(0,j,k))*dx/dy
                    + ((*Ez3D)(0,j,k+1)-(*Ez3D)(0,j,k))*dx/dz;
            }
        }
    }
    // Ex / Xmax
    if (patch->isXmax()) {
        DEBUG("Computing Xmax BC on Ex");
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ex3D)(nx_d-1,j,k) = (*Ex3D)(nx_d-2,j,k) + dx*(*rho3D)(nx_p-1,j,k)
                    - ((*Ey3D)(nx_p-1,j+1,k)-(*Ey3D)(nx_p-1,j,k))*dx/dy 
                    - ((*Ey3D)(nx_p-1,j,k+1)-(*Ez3D)(nx_p-1,j,k))*dx/dz;
            }
        }
    }

    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;

} // initE


void ElectroMagn3D::centeringE( std::vector<double> E_Add )
{
    Field3D* Ex3D  = static_cast<Field3D*>(Ex_);
    Field3D* Ey3D  = static_cast<Field3D*>(Ey_);
    Field3D* Ez3D  = static_cast<Field3D*>(Ez_);

    // Centering electrostatic fields
    for (unsigned int i=0; i<nx_d; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ex3D)(i,j,k) += E_Add[0];
            }
        }
    }
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_d; j++) {
            for (unsigned int k=0; k<nz_p; k++) {
                (*Ey3D)(i,j,k) += E_Add[1];
            }
        }
    }
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            for (unsigned int k=0; k<nz_d; k++) {
                (*Ez3D)(i,j,k) += E_Add[2];
            }
        }
    }

} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods 
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::saveMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);
    
    // Magnetic field Bx^(p,d,d)
    memcpy(&((*Bx3D_m)(0,0,0)), &((*Bx3D)(0,0,0)),nx_p*ny_d*nz_d*sizeof(double) );
    
    // Magnetic field By^(d,p,d)
    memcpy(&((*By3D_m)(0,0,0)), &((*By3D)(0,0,0)),nx_d*ny_p*nz_d*sizeof(double) );
    
    // Magnetic field Bz^(d,d,p)
    memcpy(&((*Bz3D_m)(0,0,0)), &((*Bz3D)(0,0,0)),nx_d*ny_d*nz_p*sizeof(double) );
    
}//END saveMagneticFields



//// ---------------------------------------------------------------------------------------------------------------------
//// Solve the Maxwell-Ampere equation
//// ---------------------------------------------------------------------------------------------------------------------
//void ElectroMagn3D::solveMaxwellAmpere()
//{
//    // Static-cast of the fields
//    Field3D* Ex3D = static_cast<Field3D*>(Ex_);
//    Field3D* Ey3D = static_cast<Field3D*>(Ey_);
//    Field3D* Ez3D = static_cast<Field3D*>(Ez_);
//    Field3D* Bx3D = static_cast<Field3D*>(Bx_);
//    Field3D* By3D = static_cast<Field3D*>(By_);
//    Field3D* Bz3D = static_cast<Field3D*>(Bz_);
//    Field3D* Jx3D = static_cast<Field3D*>(Jx_);
//    Field3D* Jy3D = static_cast<Field3D*>(Jy_);
//    Field3D* Jz3D = static_cast<Field3D*>(Jz_);
//    // Electric field Ex^(d,p,p)
//    for (unsigned int i=0 ; i<nx_d ; i++) {
//        for (unsigned int j=0 ; j<ny_p ; j++) {
//            for (unsigned int k=0 ; k<nz_p ; k++) {
//                (*Ex3D)(i,j,k) += -timestep*(*Jx3D)(i,j,k) + dt_ov_dy * ( (*Bz3D)(i,j+1,k) - (*Bz3D)(i,j,k) ) - dt_ov_dz * ( (*By3D)(i,j,k+1) - (*By3D)(i,j,k) );
//            }
//        }
//    }
//    
//    // Electric field Ey^(p,d,p)
//    for (unsigned int i=0 ; i<nx_p ; i++) {
//        for (unsigned int j=0 ; j<ny_d ; j++) {
//            for (unsigned int k=0 ; k<nz_p ; k++) {
//                (*Ey3D)(i,j,k) += -timestep*(*Jy3D)(i,j,k) - dt_ov_dx * ( (*Bz3D)(i+1,j,k) - (*Bz3D)(i,j,k) ) + dt_ov_dz * ( (*Bx3D)(i,j,k+1) - (*Bx3D)(i,j,k) );
//            }
//        }
//    }
//    
//    // Electric field Ez^(p,p,d)
//    for (unsigned int i=0 ;  i<nx_p ; i++) {
//        for (unsigned int j=0 ; j<ny_p ; j++) {
//            for (unsigned int k=0 ; k<nz_d ; k++) {
//                (*Ez3D)(i,j,k) += -timestep*(*Jz3D)(i,j,k) + dt_ov_dx * ( (*By3D)(i+1,j,k) - (*By3D)(i,j,k) ) - dt_ov_dy * ( (*Bx3D)(i,j+1,k) - (*Bx3D)(i,j,k) );
//            }
//        }
//    }
//
//}//END solveMaxwellAmpere


// Create a new field
Field * ElectroMagn3D::createField(string fieldname)
{
    if     (fieldname.substr(0,2)=="Ex" ) return new Field3D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Ey" ) return new Field3D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Ez" ) return new Field3D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,2)=="Bx" ) return new Field3D(dimPrim, 0, true,  fieldname);
    else if(fieldname.substr(0,2)=="By" ) return new Field3D(dimPrim, 1, true,  fieldname);
    else if(fieldname.substr(0,2)=="Bz" ) return new Field3D(dimPrim, 2, true,  fieldname);
    else if(fieldname.substr(0,2)=="Jx" ) return new Field3D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Jy" ) return new Field3D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Jz" ) return new Field3D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,3)=="Rho") return new Field3D(dimPrim, fieldname );
    
    ERROR("Cannot create field "<<fieldname);
    return NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::centerMagneticFields()
{
    // Static cast of the fields
    Field3D* Bx3D   = static_cast<Field3D*>(Bx_);
    Field3D* By3D   = static_cast<Field3D*>(By_);
    Field3D* Bz3D   = static_cast<Field3D*>(Bz_);
    Field3D* Bx3D_m = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m = static_cast<Field3D*>(Bz_m);
    
    // Magnetic field Bx^(p,d,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            for (unsigned int k=0 ; k<nz_d ; k++) {
                (*Bx3D_m)(i,j,k) = ( (*Bx3D)(i,j,k) + (*Bx3D_m)(i,j,k) )*0.5;
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            for (unsigned int k=0 ; k<nz_d ; k++) {
                (*By3D_m)(i,j,k) = ( (*By3D)(i,j,k) + (*By3D_m)(i,j,k) )*0.5;
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            for (unsigned int k=0 ; k<nz_p ; k++) {
                (*Bz3D_m)(i,j,k) = ( (*Bz3D)(i,j,k) + (*Bz3D_m)(i,j,k) )*0.5;
            }
        } // end for j
    } // end for i

    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::binomialCurrentFilter()
{
    ERROR("Binomial current filtering not yet implemented in 3D3V");
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computeTotalRhoJ()
{
    // static cast of the total currents and densities
    Field3D* Jx3D    = static_cast<Field3D*>(Jx_);
    Field3D* Jy3D    = static_cast<Field3D*>(Jy_);
    Field3D* Jz3D    = static_cast<Field3D*>(Jz_);
    Field3D* rho3D   = static_cast<Field3D*>(rho_);
    
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        if( Jx_s[ispec] ) {
            Field3D* Jx3D_s  = static_cast<Field3D*>(Jx_s[ispec]);
            for (unsigned int i=0 ; i<=nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    for (unsigned int k=0 ; k<nz_p ; k++)
                        (*Jx3D)(i,j,k) += (*Jx3D_s)(i,j,k);
        }
        if( Jy_s[ispec] ) {
            Field3D* Jy3D_s  = static_cast<Field3D*>(Jy_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<=ny_p ; j++)
                    for (unsigned int k=0 ; k<nz_p ; k++)
                        (*Jy3D)(i,j,k) += (*Jy3D_s)(i,j,k);
        }
        if( Jz_s[ispec] ) {
            Field3D* Jz3D_s  = static_cast<Field3D*>(Jz_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    for (unsigned int k=0 ; k<=nz_p ; k++)
                        (*Jz3D)(i,j,k) += (*Jz3D_s)(i,j,k);
        }
        if( rho_s[ispec] ) {
            Field3D* rho3D_s  = static_cast<Field3D*>(rho_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    for (unsigned int k=0 ; k<nz_p ; k++)
                        (*rho3D)(i,j,k) += (*rho3D_s)(i,j,k);
        }
        
    }//END loop on species ispec
//END computeTotalRhoJ
}


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computePoynting() {

    Field3D* Ex3D     = static_cast<Field3D*>(Ex_);
    Field3D* Ey3D     = static_cast<Field3D*>(Ey_);
    Field3D* Ez3D     = static_cast<Field3D*>(Ez_);
    Field3D* Bx3D_m   = static_cast<Field3D*>(Bx_m);
    Field3D* By3D_m   = static_cast<Field3D*>(By_m);
    Field3D* Bz3D_m   = static_cast<Field3D*>(Bz_m);

    if (isXmin) {
        unsigned int iEy=istart[0][Ey3D->isDual(0)];
        unsigned int iBz=istart[0][Bz3D_m->isDual(0)];
        unsigned int iEz=istart[0][Ez3D->isDual(0)];
        unsigned int iBy=istart[0][By3D_m->isDual(0)];
        
        unsigned int jEy=istart[1][Ey3D->isDual(1)];
        unsigned int jBz=istart[1][Bz3D_m->isDual(1)];
        unsigned int jEz=istart[1][Ez3D->isDual(1)];
        unsigned int jBy=istart[1][By3D_m->isDual(1)];
        
        unsigned int kEy=istart[2][Ey3D->isDual(2)];
        unsigned int kBz=istart[2][Bz3D_m->isDual(2)];
        unsigned int kEz=istart[2][Ez3D->isDual(2)];
//         unsigned int kBy=istart[2][By3D_m->isDual(2)];

        for (unsigned int j=0; j<=bufsize[1][Ez3D->isDual(1)]; j++) {
            for (unsigned int k=0; k<=bufsize[2][Ey3D->isDual(2)]; k++) {
            
                double Ey__ = 0.5 *( (*Ey3D)   (iEy, jEy+j,   kEy+k)   + (*Ey3D)  (iEy,   jEy+j+1, kEy+k) );
                double Bz__ = 0.25*( (*Bz3D_m) (iBz, jBz+j,   kBz+k)   + (*Bz3D_m)(iBz+1, jBz+j,   kBz+k)
                                     +(*Bz3D_m)(iBz, jBz+j+1, kBz+k)   + (*Bz3D_m)(iBz+1, jBz+j+1, kBz+k) );
                double Ez__ = 0.5 *( (*Ez3D)   (iEz, jEz+j  , kEz+k)   + (*Ez3D)  (iEz,   jEz+j,   kEz+k+1) );
                double By__ = 0.25*( (*By3D_m) (iBy, jBy+j,   kEz+k)   + (*By3D_m)(iBy+1, jBy+j,   kEz+k)
                                     +(*By3D_m)(iBy, jBy+j,   kEz+k+1) + (*By3D_m)(iBy+1, jBy+j,   kEz+k+1) );
            
                poynting_inst[0][0] = dy*dz*timestep*(Ey__*Bz__ - Ez__*By__);
                poynting[0][0]+= poynting_inst[0][0];

            }
        }
        
    }//if Xmin
    
    
    if (isXmax) {
        
        unsigned int iEy=istart[0][Ey3D->isDual(0)]  + bufsize[0][Ey3D->isDual(0)] -1;
        unsigned int iBz=istart[0][Bz3D_m->isDual(0)] + bufsize[0][Bz3D_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez3D->isDual(0)]  + bufsize[0][Ez3D->isDual(0)] -1;
        unsigned int iBy=istart[0][By3D_m->isDual(0)] + bufsize[0][By3D_m->isDual(0)]-1;
        
        unsigned int jEy=istart[1][Ey3D->isDual(1)];
        unsigned int jBz=istart[1][Bz3D_m->isDual(1)];
        unsigned int jEz=istart[1][Ez3D->isDual(1)];
        unsigned int jBy=istart[1][By3D_m->isDual(1)];
        
        unsigned int kEy=istart[2][Ey3D->isDual(2)];
        unsigned int kBz=istart[2][Bz3D_m->isDual(2)];
        unsigned int kEz=istart[2][Ez3D->isDual(2)];
//         unsigned int kBy=istart[2][By3D_m->isDual(2)];
        
        for (unsigned int j=0; j<=bufsize[1][Ez3D->isDual(1)]; j++) {
            for (unsigned int k=0; k<=bufsize[2][Ey3D->isDual(2)]; k++) {
            
                double Ey__ = 0.5 *( (*Ey3D)   (iEy, jEy+j,   kEy+k)   + (*Ey3D)  (iEy,   jEy+j+1, kEy+k) );
                double Bz__ = 0.25*( (*Bz3D_m) (iBz, jBz+j,   kBz+k)   + (*Bz3D_m)(iBz+1, jBz+j,   kBz+k)
                                     +(*Bz3D_m)(iBz, jBz+j+1, kBz+k)   + (*Bz3D_m)(iBz+1, jBz+j+1, kBz+k) );
                double Ez__ = 0.5 *( (*Ez3D)   (iEz, jEz+j,   kEz+k)   + (*Ez3D)  (iEz,   jEz+j,   kEz+k+1));
                double By__ = 0.25*( (*By3D_m) (iBy, jBy+j,   kEz+k)   + (*By3D_m)(iBy+1, jBy+j,   kEz+k)
                                     +(*By3D_m)(iBy, jBy+j,   kEz+k+1) + (*By3D_m)(iBy+1, jBy+j,   kEz+k+1) );
            
                poynting_inst[1][0] = dy*dz*timestep*(Ey__*Bz__ - Ez__*By__);
                poynting[1][0]+= poynting_inst[1][0];

            }
        }
        
    }//if Xmax
    
    if (isYmin) {
        
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)]; 
        
        unsigned int jEz=istart[1][Ez_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBz=istart[1][Bz_m->isDual(1)];

        unsigned int kEz=istart[2][Ez_->isDual(2)];
        unsigned int kBx=istart[2][Bx_m->isDual(2)];
        unsigned int kEx=istart[2][Ex_->isDual(2)];
//         unsigned int kBz=istart[2][Bz_m->isDual(2)];
        
        for (unsigned int i=0; i<=bufsize[0][Ez3D->isDual(0)]; i++) {
            for (unsigned int k=0; k<=bufsize[2][Ex3D->isDual(2)]; k++) {
                double Ez__ = 0.5 *( (*Ez3D)   (iEz+i, jEz,   kEz+k)   + (*Ez3D)  (iEz+i,   jEz,   kEz+k+1) );
                double Bx__ = 0.25*( (*Bx3D_m) (iBx+i, jBx,   kBx+k)   + (*Bx3D_m)(iBx+i,   jBx+1, kBx+k) 
                                     +(*Bx3D_m)(iBx+i, jBx,   kBx+k+1) + (*Bx3D_m)(iBx+i,   jBx+1, kBx+k+1) );
                double Ex__ = 0.5 *( (*Ex3D)   (iEx+i, jEx,   kEx+k)   + (*Ex3D)  (iEx+i+1, jEx,   kEx+k) );
                double Bz__ = 0.25*( (*Bz3D_m) (iBz+i, jBz,   kEx+k)   + (*Bz3D_m)(iBz+i+1, jBz,   kEx+k)
                                     +(*Bz3D_m)(iBz+i, jBz+1, kEx+k)   + (*Bz3D_m)(iBz+i+1, jBz+1, kEx+k) );
            
                poynting_inst[0][1] = dx*dz*timestep*(Ez__*Bx__ - Ex__*Bz__);
                poynting[0][1] += poynting_inst[0][1];
            }
        }

    }// if Ymin
    
    if (isYmax) {

        unsigned int iEz=istart[0][Ez3D->isDual(0)];
        unsigned int iBx=istart[0][Bx3D_m->isDual(0)];
        unsigned int iEx=istart[0][Ex3D->isDual(0)];
        unsigned int iBz=istart[0][Bz3D_m->isDual(0)];
        
        unsigned int jEz=istart[1][Ez3D->isDual(1)]  + bufsize[1][Ez3D->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx3D_m->isDual(1)] + bufsize[1][Bx3D_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex3D->isDual(1)]  + bufsize[1][Ex3D->isDual(1)] -1;
        unsigned int jBz=istart[1][Bz3D_m->isDual(1)] + bufsize[1][Bz3D_m->isDual(1)]-1;
        
        unsigned int kEz=istart[2][Ez_->isDual(2)];
        unsigned int kBx=istart[2][Bx_m->isDual(2)];
        unsigned int kEx=istart[2][Ex_->isDual(2)];
//         unsigned int kBz=istart[2][Bz_m->isDual(2)];
        
        for (unsigned int i=0; i<=bufsize[0][Ez3D->isDual(0)]; i++) {
            for (unsigned int k=0; k<=bufsize[2][Ex3D->isDual(2)]; k++) {
                double Ez__ = 0.5 *( (*Ez3D)   (iEz+i, jEz,   kEz+k)   + (*Ez3D)  (iEz+i,   jEz,   kEz+k+1) );
                double Bx__ = 0.25*( (*Bx3D_m) (iBx+i, jBx,   kBx+k)   + (*Bx3D_m)(iBx+i,   jBx+1, kBx+k) 
                                     +(*Bx3D_m)(iBx+i, jBx,   kBx+k+1) + (*Bx3D_m)(iBx+i,   jBx+1, kBx+k+1) );
                double Ex__ = 0.5 *( (*Ex3D)   (iEx+i, jEx,   kEx+k)   + (*Ex3D)  (iEx+i+1, jEx,   kEx+k) );
                double Bz__ = 0.25*( (*Bz3D_m) (iBz+i, jBz,   kEx+k)   + (*Bz3D_m)(iBz+i+1, jBz,   kEx+k)
                                     +(*Bz3D_m)(iBz+i, jBz+1, kEx+k)   + (*Bz3D_m)(iBz+i+1, jBz+1, kEx+k) );
            
                poynting_inst[1][1] = dx*dz*timestep*(Ez__*Bx__ - Ex__*Bz__);
                poynting[1][1] += poynting_inst[1][1];
            }
        }

    }//if Ymax

    if (isZmin) {
        
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)]; 
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];
        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];

        unsigned int kEx=istart[2][Ex_->isDual(2)];
        unsigned int kBy=istart[2][By_m->isDual(2)];
        unsigned int kEy=istart[2][Ey_->isDual(2)];
//         unsigned int kBx=istart[2][Bx_m->isDual(2)];

        for (unsigned int i=0; i<=bufsize[0][Ez3D->isDual(0)]; i++) {
            for (unsigned int j=0; j<=bufsize[1][Ex3D->isDual(1)]; j++) {

                double Ex__ = 0.5 *( (*Ex3D)   (iEx+i, jEx+j, kEx)   + (*Ex3D)  (iEx+i+1, jEx+j,   kEx) );
                double By__ = 0.25*( (*By3D_m) (iBy+i, jBy+j, kBy)   + (*By3D_m)(iBy+i+1, jBy+j,   kBy) 
                                     +(*By3D_m)(iBy+i, jBy+j, kBy+1) + (*By3D_m)(iBy+i+1, jBy+j,   kBy+1) );
                double Ey__ = 0.5 *( (*Ey3D)   (iEy+i, jEy+j, kEy)   + (*Ey3D)  (iEy+i,   jEy+j+1, kEy) );
                double Bx__ = 0.25*( (*Bx3D_m) (iBx+i, jBx+j, kEx)   + (*Bx3D_m)(iBx+i,   jBx+j+1, kEx)
                                     +(*Bx3D_m)(iBx+i, jBx+j, kEx+1) + (*Bx3D_m)(iBx+i,   jBx+j+1, kEx+1) );
            
                poynting_inst[0][2] = dx*dy*timestep*(Ex__*By__ - Ey__*Bx__);
                poynting[0][2] += poynting_inst[0][2];
            }
        }

    }

    if (isZmax) {
        
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)]; 
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBy=istart[1][By_m->isDual(1)];
        unsigned int jEy=istart[1][Ey_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];

        unsigned int kEx=istart[2][Ex_->isDual(2)] + bufsize[2][Ex_->isDual(2)] - 1;
        unsigned int kBy=istart[2][By_m->isDual(2)] + bufsize[2][By_m->isDual(2)] - 1;
        unsigned int kEy=istart[2][Ey_->isDual(2)] + bufsize[2][Ey_->isDual(2)] - 1;
//         unsigned int kBx=istart[2][Bx_m->isDual(2)] + bufsize[2][Bx_m->isDual(2)] - 1;

        for (unsigned int i=0; i<=bufsize[0][Ez3D->isDual(0)]; i++) {
            for (unsigned int j=0; j<=bufsize[1][Ex3D->isDual(1)]; j++) {

                double Ex__ = 0.5 *( (*Ex3D)   (iEx+i, jEx+j, kEx)   + (*Ex3D)  (iEx+i+1, jEx+j,   kEx) );
                double By__ = 0.25*( (*By3D_m) (iBy+i, jBy+j, kBy)   + (*By3D_m)(iBy+i+1, jBy+j,   kBy) 
                                     +(*By3D_m)(iBy+i, jBy+j, kBy+1) + (*By3D_m)(iBy+i+1, jBy+j,   kBy+1) );
                double Ey__ = 0.5 *( (*Ey3D)   (iEy+i, jEy+j, kEy)   + (*Ey3D)  (iEy+i,   jEy+j+1, kEy) );
                double Bx__ = 0.25*( (*Bx3D_m) (iBx+i, jBx+j, kEx)   + (*Bx3D_m)(iBx+i,   jBx+j+1, kEx)
                                     +(*Bx3D_m)(iBx+i, jBx+j, kEx+1) + (*Bx3D_m)(iBx+i,   jBx+j+1, kEx+1) );
            
                poynting_inst[1][2] = dx*dy*timestep*(Ex__*By__ - Ey__*Bx__);
                poynting[1][2] += poynting_inst[1][2];
            }
        }

    }

}

void ElectroMagn3D::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    
    Field3D* field3D=static_cast<Field3D*>(my_field);
    
    vector<double> pos(3);
    pos[0]      = dx*((double)(patch->getCellStartingGlobalIndex(0))+(field3D->isDual(0)?-0.5:0.));
    double pos1 = dy*((double)(patch->getCellStartingGlobalIndex(1))+(field3D->isDual(1)?-0.5:0.));
    double pos2 = dz*((double)(patch->getCellStartingGlobalIndex(2))+(field3D->isDual(2)?-0.5:0.));
    int N0 = (int)field3D->dims()[0];
    int N1 = (int)field3D->dims()[1];
    int N2 = (int)field3D->dims()[2];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (int i=0 ; i<N0 ; i++) {
        pos[1] = pos1;
        for (int j=0 ; j<N1 ; j++) {
            pos[2] = pos2;
            for (int k=0 ; k<N2 ; k++) {
                (*field3D)(i,j,k) += profile->valueAt(pos);
                pos[2] += dz;
            }
            pos[1] += dy;
        }
        pos[0] += dx;
    }
    
    for (auto& embc: emBoundCond) {
        if (embc) embc->save_fields(my_field, patch);
    }
}



void ElectroMagn3D::initAntennas(Patch* patch)
{
    
    // Filling the space profiles of antennas
    for (unsigned int i=0; i<antennas.size(); i++) {
        if      (antennas[i].fieldName == "Jx")
            antennas[i].field = new Field3D(dimPrim, 0, false, "Jx");
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = new Field3D(dimPrim, 1, false, "Jy");
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = new Field3D(dimPrim, 2, false, "Jz");
        else {
            ERROR("Antenna cannot be applied to field "<<antennas[i].fieldName);
        }
        
        if (antennas[i].field) 
            applyExternalField(antennas[i].field, antennas[i].space_profile, patch);
    }

}

