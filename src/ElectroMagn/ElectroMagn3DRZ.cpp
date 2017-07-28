#include "ElectroMagn3DRZ.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "Params.h"
#include "Field2D.h"
#include "cField2D.h"

#include "Patch.h"
#include <cstring>

#include "Profile.h"

#include "ElectroMagnBC.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn3DRZ
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3DRZ::ElectroMagn3DRZ(Params &params, vector<Species*>& vecSpecies, Patch* patch) : 
  ElectroMagn(params, vecSpecies, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax())
{    
    
    initElectroMagn3DRZQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field2D(("Jx_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        Jy_s[ispec]  = new Field2D(("Jy_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        Jz_s[ispec]  = new Field2D(("Jz_" +vecSpecies[ispec]->species_type).c_str(), dimPrim);
        rho_s[ispec] = new Field2D(("Rho_"+vecSpecies[ispec]->species_type).c_str(), dimPrim);
    }
    
}//END constructor Electromagn3D


ElectroMagn3DRZ::ElectroMagn3DRZ( ElectroMagn3DRZ* emFields, Params &params, Patch* patch ) : 
    ElectroMagn(emFields, params, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax())
{
    
    initElectroMagn3DRZQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int ispec=0; ispec<n_species*2*nmodes; ispec++) {              // *2 for real and imaginary components, *nmodes for number of modes
        if ( emFields->Jx_s[ispec] != NULL ) {
            if ( emFields->Jx_s[ispec]->data_ != NULL )
                Jx_s[ispec]  = new Field2D(dimPrim, 0, false, emFields->Jx_s[ispec]->name);
            else
                Jx_s[ispec]  = new Field2D(emFields->Jx_s[ispec]->name, dimPrim);
        }
        if ( emFields->Jy_s[ispec] != NULL ) {
            if ( emFields->Jy_s[ispec]->data_ != NULL )
                Jy_s[ispec]  = new Field2D(dimPrim, 1, false, emFields->Jy_s[ispec]->name);
            else
                Jy_s[ispec]  = new Field2D(emFields->Jy_s[ispec]->name, dimPrim);
        }
        if ( emFields->Jz_s[ispec] != NULL ) {
            if ( emFields->Jz_s[ispec]->data_ != NULL )
                Jz_s[ispec]  = new Field2D(dimPrim, 2, false, emFields->Jz_s[ispec]->name);
            else
                Jz_s[ispec]  = new Field2D(emFields->Jz_s[ispec]->name, dimPrim);
        }
        if ( emFields->rho_s[ispec] != NULL ) {
            if ( emFields->rho_s[ispec]->data_ != NULL )
                rho_s[ispec] = new Field2D(dimPrim, emFields->rho_s[ispec]->name );
            else
                rho_s[ispec]  = new Field2D(emFields->rho_s[ispec]->name, dimPrim);
        }
    }


}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn3D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::initElectroMagn3DRZQuantities(Params &params, Patch* patch)
{

    // Species charge currents and density
    Jx_s.resize(n_species*2*nmodes);
    Jy_s.resize(n_species*2*nmodes);
    Jz_s.resize(n_species*2*nmodes);
    rho_s.resize(n_species*2*nmodes);
    for (unsigned int ispec=0; ispec<n_species*2*nmodes; ispec++) {
        Jx_s[ispec]  = NULL;
        Jy_s[ispec]  = NULL;
        Jz_s[ispec]  = NULL;
        rho_s[ispec] = NULL;
    }

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
// Destructor for Electromagn3DRZ
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3DRZ::~ElectroMagn3DRZ()
{
}//END ElectroMagn3DRZ



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


void ElectroMagn3DRZ::initPoisson(Patch *patch)
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
    if (patch->isXmin()) {
        index_min_p_[0] = 0;
    }
    if (patch->isXmax()) {
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

double ElectroMagn3DRZ::compute_r()
{
    double rnew_dot_rnew_local(0.);
    for (unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++) {
        for (unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++) {
            rnew_dot_rnew_local += (*r_)(i,j)*(*r_)(i,j);
        }
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn3DRZ::compute_Ap(Patch* patch)
{
} // compute_pAp

double ElectroMagn3DRZ::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagn3DRZ::update_pand_r(double r_dot_r, double p_dot_Ap)
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for(unsigned int i=0; i<nx_p; i++) {
        for(unsigned int j=0; j<ny_p; j++) {
            (*phi_)(i,j) += alpha_k * (*p_)(i,j);
            (*r_)(i,j)   -= alpha_k * (*Ap_)(i,j);
        }
    }

} // update_pand_r

void ElectroMagn3DRZ::update_p(double rnew_dot_rnew, double r_dot_r)
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for (unsigned int i=0; i<nx_p; i++) {
        for(unsigned int j=0; j<ny_p; j++) {
            (*p_)(i,j) = (*r_)(i,j) + beta_k * (*p_)(i,j);
        }
    }
} // update_p

void ElectroMagn3DRZ::initE(Patch *patch)
{

    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;

} // initE


void ElectroMagn3DRZ::centeringE( std::vector<double> E_Add )
{
    Field2D* Ex2D  = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D  = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D  = static_cast<Field2D*>(Ez_);

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
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ez2D)(i,j) += E_Add[2];
        }
    }

} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods 
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::saveMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d,d)
    memcpy(&((*Bx2D_m)(0,0)), &((*Bx2D)(0,0)),nx_p*ny_d*sizeof(double) );
    
    // Magnetic field By^(d,p,d)
    memcpy(&((*By2D_m)(0,0)), &((*By2D)(0,0)),nx_d*ny_p*sizeof(double) );
    
    // Magnetic field Bz^(d,d,p)
    memcpy(&((*Bz2D_m)(0,0)), &((*Bz2D)(0,0)),nx_d*ny_d*sizeof(double) );

}//END saveMagneticFields



//// ---------------------------------------------------------------------------------------------------------------------
//// Solve the Maxwell-Ampere equation
//// ---------------------------------------------------------------------------------------------------------------------
//void ElectroMagn3DRZ::solveMaxwellAmpere()
//{
//    // Static-cast of the fields
//    Field2D* Ex2D = static_cast<Field2D*>(Ex_);
//    Field2D* Ey2D = static_cast<Field2D*>(Ey_);
//    Field2D* Ez2D = static_cast<Field2D*>(Ez_);
//    Field2D* Bx2D = static_cast<Field2D*>(Bx_);
//    Field2D* By2D = static_cast<Field2D*>(By_);
//    Field2D* Bz2D = static_cast<Field2D*>(Bz_);
//    Field2D* Jx2D = static_cast<Field2D*>(Jx_);
//    Field2D* Jy2D = static_cast<Field2D*>(Jy_);
//    Field2D* Jz2D = static_cast<Field2D*>(Jz_);
//    // Electric field Ex^(d,p,p)
//    for (unsigned int i=0 ; i<nx_d ; i++) {
//        for (unsigned int j=0 ; j<ny_p ; j++) {
//            for (unsigned int k=0 ; k<nz_p ; k++) {
//                (*Ex2D)(i,j) += -timestep*(*Jx2D)(i,j) + dt_ov_dy * ( (*Bz2D)(i,j+1) - (*Bz2D)(i,j) ) - dt_ov_dz * ( (*By2D)(i,j) - (*By2D)(i,j) );
//            }
//        }
//    }
//    
//    // Electric field Ey^(p,d,p)
//    for (unsigned int i=0 ; i<nx_p ; i++) {
//        for (unsigned int j=0 ; j<ny_d ; j++) {
//            for (unsigned int k=0 ; k<nz_p ; k++) {
//                (*Ey2D)(i,j) += -timestep*(*Jy2D)(i,j) - dt_ov_dx * ( (*Bz2D)(i+1,j) - (*Bz2D)(i,j) ) + dt_ov_dz * ( (*Bx2D)(i,j) - (*Bx2D)(i,j) );
//            }
//        }
//    }
//    
//    // Electric field Ez^(p,p,d)
//    for (unsigned int i=0 ;  i<nx_p ; i++) {
//        for (unsigned int j=0 ; j<ny_p ; j++) {
//            for (unsigned int k=0 ; k<nz_d ; k++) {
//                (*Ez2D)(i,j) += -timestep*(*Jz2D)(i,j) + dt_ov_dx * ( (*By2D)(i+1,j) - (*By2D)(i,j) ) - dt_ov_dy * ( (*Bx2D)(i,j+1) - (*Bx2D)(i,j) );
//            }
//        }
//    }
//
//}//END solveMaxwellAmpere


// Create a new field
Field * ElectroMagn3DRZ::createField(string fieldname)
{
    if     (fieldname.substr(0,2)=="Ex" ) return new Field2D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Ey" ) return new Field2D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Ez" ) return new Field2D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,2)=="Bx" ) return new Field2D(dimPrim, 0, true,  fieldname);
    else if(fieldname.substr(0,2)=="By" ) return new Field2D(dimPrim, 1, true,  fieldname);
    else if(fieldname.substr(0,2)=="Bz" ) return new Field2D(dimPrim, 2, true,  fieldname);
    else if(fieldname.substr(0,2)=="Jx" ) return new Field2D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Jy" ) return new Field2D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Jz" ) return new Field2D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,3)=="Rho") return new Field2D(dimPrim, fieldname );
    
    ERROR("Cannot create field "<<fieldname);
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::centerMagneticFields()
{
    // Static cast of the fields
    Field2D* Bx2D   = static_cast<Field2D*>(Bx_);
    Field2D* By2D   = static_cast<Field2D*>(By_);
    Field2D* Bz2D   = static_cast<Field2D*>(Bz_);
    Field2D* Bx2D_m = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m = static_cast<Field2D*>(Bz_m);
    
    // Magnetic field Bx^(p,d,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bx2D_m)(i,j) = ( (*Bx2D)(i,j) + (*Bx2D_m)(i,j) )*0.5;
        }
    }
    
    // Magnetic field By^(d,p,d)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*By2D_m)(i,j) = ( (*By2D)(i,j) + (*By2D_m)(i,j) )*0.5;
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*Bz2D_m)(i,j) = ( (*Bz2D)(i,j) + (*Bz2D_m)(i,j) )*0.5;
        } // end for j
    } // end for i

    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::binomialCurrentFilter()
{
    ERROR("Binomial current filtering not yet implemented in 2D3V");
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computeTotalRhoJ()
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
        if( Jx_s[ispec] ) {
            Field2D* Jx2D_s  = static_cast<Field2D*>(Jx_s[ispec]);
            for (unsigned int i=0 ; i<=nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    (*Jx2D)(i,j) += (*Jx2D_s)(i,j);
        }
        if( Jy_s[ispec] ) {
            Field2D* Jy2D_s  = static_cast<Field2D*>(Jy_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<=ny_p ; j++)
                    (*Jy2D)(i,j) += (*Jy2D_s)(i,j);
        }
        if( Jz_s[ispec] ) {
            Field2D* Jz2D_s  = static_cast<Field2D*>(Jz_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    (*Jz2D)(i,j) += (*Jz2D_s)(i,j);
        }
        if( rho_s[ispec] ) {
            Field2D* rho2D_s  = static_cast<Field2D*>(rho_s[ispec]);
            for (unsigned int i=0 ; i<nx_p ; i++)
                for (unsigned int j=0 ; j<ny_p ; j++)
                    (*rho2D)(i,j) += (*rho2D_s)(i,j);
        }
        
    }//END loop on species ispec
//END computeTotalRhoJ
}


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computePoynting() {

    Field2D* Ex2D     = static_cast<Field2D*>(Ex_);
    Field2D* Ey2D     = static_cast<Field2D*>(Ey_);
    Field2D* Ez2D     = static_cast<Field2D*>(Ez_);
    Field2D* Bx2D_m   = static_cast<Field2D*>(Bx_m);
    Field2D* By2D_m   = static_cast<Field2D*>(By_m);
    Field2D* Bz2D_m   = static_cast<Field2D*>(Bz_m);

    if (isXmin) {
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
        
    }//if Xmin
    
    
    if (isXmax) {
        
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
            poynting[1][0]+= poynting_inst[1][0];

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

        for (unsigned int i=0; i<=bufsize[0][Ez2D->isDual(0)]; i++) {
            double Ez__ = (*Ez2D)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
            
            poynting_inst[0][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[0][1] += poynting_inst[0][1];
        }

    }// if Ymin
    
    if (isYmax) {

        unsigned int iEz=istart[0][Ez2D->isDual(0)];
        unsigned int iBx=istart[0][Bx2D_m->isDual(0)];
        unsigned int iEx=istart[0][Ex2D->isDual(0)];
        unsigned int iBz=istart[0][Bz2D_m->isDual(0)];
        
        unsigned int jEz=istart[1][Ez2D->isDual(1)]  + bufsize[1][Ez2D->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx2D_m->isDual(1)] + bufsize[1][Bx2D_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex2D->isDual(1)]  + bufsize[1][Ex2D->isDual(1)] -1;
        unsigned int jBz=istart[1][Bz2D_m->isDual(1)] + bufsize[1][Bz2D_m->isDual(1)]-1;
        
        for (unsigned int i=0; i<=bufsize[0][Ez2D->isDual(0)]; i++) {
            double Ez__ = (*Ez2D)(iEz+i,jEz);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
            
            poynting_inst[1][1] = dx*timestep*(Ez__*Bx__ - Ex__*Bz__);
            poynting[1][1] += poynting_inst[1][1];
        }

    }//if Ymax

}

void ElectroMagn3DRZ::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    
    Field2D* field2D=static_cast<Field2D*>(my_field);
    
    vector<double> pos(2);
    pos[0]      = dx*((double)(patch->getCellStartingGlobalIndex(0))+(field2D->isDual(0)?-0.5:0.));
    double pos1 = dy*((double)(patch->getCellStartingGlobalIndex(1))+(field2D->isDual(1)?-0.5:0.));
    int N0 = (int)field2D->dims()[0];
    int N1 = (int)field2D->dims()[1];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (int i=0 ; i<N0 ; i++) {
        pos[1] = pos1;
        for (int j=0 ; j<N1 ; j++) {
            (*field2D)(i,j) += profile->valueAt(pos);
            pos[1] += dy;
        }
        pos[0] += dx;
    }
    
    for (auto& embc: emBoundCond) {
        if (embc) embc->save_fields(my_field, patch);
    }
}



void ElectroMagn3DRZ::initAntennas(Patch* patch)
{
    
    // Filling the space profiles of antennas
    for (unsigned int i=0; i<antennas.size(); i++) {
        if      (antennas[i].fieldName == "Jx")
            antennas[i].field = new Field2D(dimPrim, 0, false, "Jx");
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = new Field2D(dimPrim, 1, false, "Jy");
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = new Field2D(dimPrim, 2, false, "Jz");
        else {
            ERROR("Antenna cannot be applied to field "<<antennas[i].fieldName);
        }
        
        if (antennas[i].field) 
            applyExternalField(antennas[i].field, antennas[i].space_profile, patch);
    }

}

