#include "ElectroMagn1D.h"

#include <cmath>

#include <sstream>
#include <string>
#include <iostream>

#include "Params.h"
#include "Field1D.h"

#include "Patch.h"

#include "Profile.h"
#include "MF_Solver1D_Yee.h"

#include "ElectroMagnBC.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::ElectroMagn1D(Params &params, vector<Species*>& vecSpecies, Patch* patch)
  : ElectroMagn(params, vecSpecies, patch),
isWestern(patch->isWestern()),
isEastern(patch->isEastern())
{
    oversize_ = oversize[0];
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step
    dx       = cell_length[0];
    dt_ov_dx = timestep/cell_length[0];
    dx_ov_dt = 1.0/dt_ov_dx;
    
    // Electromagnetic fields
    // ----------------------
    // number of nodes of the primal-grid
    nx_p = n_space[0]+1 + 2*oversize[0];
    // number of nodes of the dual-grid
    nx_d = n_space[0]+2 + 2*oversize[0];
    // dimPrim/dimDual = nx_p/nx_d
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    for (size_t i=0 ; i<nDim_field ; i++) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    
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
    
    // for (unsigned int i=0 ; i<nx_d ; i++) {
    //         double x = ( (double)(smpi1D->getCellStartingGlobalIndex(0)+i-0.5) )*params.cell_length[0];
    //         (*By_)(i) = 0.001 * sin(x * 2.0*M_PI/params.sim_length[0] * 40.0);
    //     }
    //     smpi1D->exchangeField(By_);
    //     for (unsigned int i=0 ; i<nx_d ; i++) {
    // //        double x = ( (double)(smpi1D->getCellStartingGlobalIndex(0)+i-0.5) )*params.cell_length[0];
    //         (*By_m)(i) = (*By_)(i);
    //     }
    //     
    // Allocation of time-averaged EM fields
    Ex_avg  = new Field1D(dimPrim, 0, false, "Ex_avg");
    Ey_avg  = new Field1D(dimPrim, 1, false, "Ey_avg");
    Ez_avg  = new Field1D(dimPrim, 2, false, "Ez_avg");
    Bx_avg  = new Field1D(dimPrim, 0, true,  "Bx_avg");
    By_avg  = new Field1D(dimPrim, 1, true,  "By_avg");
    Bz_avg  = new Field1D(dimPrim, 2, true,  "Bz_avg");
    
    // Total charge currents and densities
    Jx_   = new Field1D(dimPrim, 0, false, "Jx");
    Jy_   = new Field1D(dimPrim, 1, false, "Jy");
    Jz_   = new Field1D(dimPrim, 2, false, "Jz");
    rho_  = new Field1D(dimPrim, "Rho" );
    
    // Charge currents currents and density for each species
    
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Jx_s[ispec]  = new Field1D(dimPrim, 0, false, ("Jx_"+vecSpecies[ispec]->species_type).c_str());
        Jy_s[ispec]  = new Field1D(dimPrim, 1, false, ("Jy_"+vecSpecies[ispec]->species_type).c_str());
        Jz_s[ispec]  = new Field1D(dimPrim, 2, false, ("Jz_"+vecSpecies[ispec]->species_type).c_str());
        rho_s[ispec] = new Field1D(dimPrim, ("Rho_"+vecSpecies[ispec]->species_type).c_str());
    }
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for (size_t i=0 ; i<nDim_field ; i++) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    
    
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
                    if ( (patch->Pcoordinates[i]!=0) && (patch->Pcoordinates[i]!=(unsigned int)params.number_of_patches[i]-1) ) 
                        bufsize[i][isDual]--;
                }
                
            } // if ( params.number_of_patches[i]!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field 
    
}//END constructor Electromagn1D



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn1D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn1D::~ElectroMagn1D()
{
}

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

void ElectroMagn1D::initPoisson(Patch *patch)
{
    Field1D* rho1D = static_cast<Field1D*>(rho_);
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------
    
    index_min_p_.resize(1,0);
    index_max_p_.resize(1,0);
    
    index_min_p_[0] = oversize[0];
    index_max_p_[0] = nx_p - 2 - oversize[0];
    if (patch->isWestern()) {
        index_min_p_[0] = 0;
    }
    if (patch->isEastern()) {
        index_max_p_[0] = nx_p-1;
    }
    
    phi_ = new Field1D(dimPrim);    // scalar potential
    r_   = new Field1D(dimPrim);    // residual vector
    p_   = new Field1D(dimPrim);    // direction vector
    Ap_  = new Field1D(dimPrim);    // A*p vector
    
    double       dx_sq          = dx*dx;
    
    // phi: scalar potential, r: residual and p: direction
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*phi_)(i)   = 0.0;
        (*r_)(i)     = -dx_sq * (*rho1D)(i);
        (*p_)(i)     = (*r_)(i);
    }
} // initPoisson

double ElectroMagn1D::compute_r()
{
    double rnew_dot_rnew_local(0.);
    for (unsigned int i=index_min_p_[0] ; i<=index_max_p_[0] ; i++)
        rnew_dot_rnew_local += (*r_)(i)*(*r_)(i);
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn1D::compute_Ap(Patch *patch)
{
    // vector product Ap = A*p
    for (unsigned int i=1 ; i<dimPrim[0]-1 ; i++)
        (*Ap_)(i) = (*p_)(i-1) - 2.0*(*p_)(i) + (*p_)(i+1);
        
    // apply BC on Ap
    if (patch->isWestern()) (*Ap_)(0)      = (*p_)(1)      - 2.0*(*p_)(0);
    if (patch->isEastern()) (*Ap_)(nx_p-1) = (*p_)(nx_p-2) - 2.0*(*p_)(nx_p-1); 
    
} // compute_Ap

double ElectroMagn1D::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    for (unsigned int i=index_min_p_[0] ; i<=index_max_p_[0] ; i++)
        p_dot_Ap_local += (*p_)(i)*(*Ap_)(i);
    return p_dot_Ap_local;

} // compute_pAp

void ElectroMagn1D::update_pand_r(double r_dot_r, double p_dot_Ap)
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*phi_)(i) += alpha_k * (*p_)(i);
        (*r_)(i)   -= alpha_k * (*Ap_)(i);
    }

} // update_pand_r

void ElectroMagn1D::update_p(double rnew_dot_rnew, double r_dot_r)
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for (unsigned int i=0 ; i<dimPrim[0] ; i++)
        (*p_)(i) = (*r_)(i) + beta_k * (*p_)(i);
} // update_p

void ElectroMagn1D::initE(Patch *patch)
{
    Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
    Field1D* rho1D = static_cast<Field1D*>(rho_);

    // ----------------------------------
    // Compute the electrostatic field Ex
    // ----------------------------------
    
    for (unsigned int i=1; i<nx_d-1; i++)
        (*Ex1D)(i) = ((*phi_)(i-1)-(*phi_)(i))/dx;
    
    // BC on Ex
    if (patch->isWestern()) (*Ex1D)(0)      = (*Ex1D)(1)      - dx*(*rho1D)(0);
    if (patch->isEastern()) (*Ex1D)(nx_d-1) = (*Ex1D)(nx_d-2) + dx*(*rho1D)(nx_p-1);
    
    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;

} // initE

void ElectroMagn1D::centeringE( std::vector<double> E_Add )
{
    Field1D* Ex1D  = static_cast<Field1D*>(Ex_);
    for (unsigned int i=0; i<nx_d; i++)
        (*Ex1D)(i) += E_Add[0];

} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods 
// ---------------------------------------------------------------------------------------------------------------------


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
        (*Ex1D)(ix)= (*Ex1D)(ix) - timestep* (*Jx1D)(ix) ;
    }
    // Transverse fields ey, ez  are defined on the primal grid
    //for (unsigned int ix=0 ; ix<nx_p ; ix++) {
    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
        (*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - timestep * (*Jy1D)(ix) ;
        (*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - timestep * (*Jz1D)(ix) ;
    }
    
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
// Reset/Increment the averaged fields
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::incrementAvgFields(unsigned int time_step)
{
    // Static cast of the fields
    Field1D* Ex1D     = static_cast<Field1D*>(Ex_);
    Field1D* Ey1D     = static_cast<Field1D*>(Ey_);
    Field1D* Ez1D     = static_cast<Field1D*>(Ez_);
    Field1D* Bx1D_m   = static_cast<Field1D*>(Bx_m);
    Field1D* By1D_m   = static_cast<Field1D*>(By_m);
    Field1D* Bz1D_m   = static_cast<Field1D*>(Bz_m);
    Field1D* Ex1D_avg = static_cast<Field1D*>(Ex_avg);
    Field1D* Ey1D_avg = static_cast<Field1D*>(Ey_avg);
    Field1D* Ez1D_avg = static_cast<Field1D*>(Ez_avg);
    Field1D* Bx1D_avg = static_cast<Field1D*>(Bx_avg);
    Field1D* By1D_avg = static_cast<Field1D*>(By_avg);
    Field1D* Bz1D_avg = static_cast<Field1D*>(Bz_avg);
    
    // for Ey^(p), Ez^(p) & Bx^(p)
    for (unsigned int i=0 ; i<dimPrim[0] ; i++) {
        (*Ey1D_avg)(i) += (*Ey1D)(i);
        (*Ez1D_avg)(i) += (*Ez1D)(i);
        (*Bx1D_avg)(i) += (*Bx1D_m)(i);
    }
    
    // for Ex^(d), By^(d) & Bz^(d)
    for (unsigned int i=0 ; i<dimDual[0] ; i++) {
        (*Ex1D_avg)(i) += (*Ex1D)(i);
        (*By1D_avg)(i) += (*By1D_m)(i);
        (*Bz1D_avg)(i) += (*Bz1D_m)(i);
    }
    
}//END incrementAvgFields



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalRhoJ()
{
    Field1D* Jx1D    = static_cast<Field1D*>(Jx_);
    Field1D* Jy1D    = static_cast<Field1D*>(Jy_);
    Field1D* Jz1D    = static_cast<Field1D*>(Jz_);
    Field1D* rho1D   = static_cast<Field1D*>(rho_);
    
    for (unsigned int ispec=0; ispec<n_species; ispec++) {
        Field1D* Jx1D_s  = static_cast<Field1D*>(Jx_s[ispec]);
        Field1D* Jy1D_s  = static_cast<Field1D*>(Jy_s[ispec]);
        Field1D* Jz1D_s  = static_cast<Field1D*>(Jz_s[ispec]);
        Field1D* rho1D_s = static_cast<Field1D*>(rho_s[ispec]);
        
        for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
            (*Jx1D)(ix)  += (*Jx1D_s)(ix);
            (*Jy1D)(ix)  += (*Jy1D_s)(ix);
            (*Jz1D)(ix)  += (*Jz1D_s)(ix);
            (*rho1D)(ix) += (*rho1D_s)(ix);
        }

        {
            (*Jx1D)(dimPrim[0])  += (*Jx1D_s)(dimPrim[0]);
        }
    }//END loop on species ispec
}

// --------------------------------------------------------------------------
// Compute Poynting (return the electromagnetic energy injected at the border
// --------------------------------------------------------------------------
void ElectroMagn1D::computePoynting() {
    
    // Western border (Energy injected = +Poynting)
    if (isWestern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)];
        unsigned int iBz=istart[0][Bz_m->isDual(0)];
        unsigned int iEz=istart[0][Ez_->isDual(0)];
        unsigned int iBy=istart[0][By_m->isDual(0)];
        
        poynting_inst[0][0]=0.5*timestep*((*Ey_)(iEy) * ((*Bz_m)(iBz) + (*Bz_m)(iBz+1)) -
                                          (*Ez_)(iEz) * ((*By_m)(iBy) + (*By_m)(iBy+1)));
        poynting[0][0] += poynting_inst[0][0];
    }
    
    // Eastern border (Energy injected = -Poynting)
    if (isEastern) {
        unsigned int iEy=istart[0][Ey_->isDual(0)]  + bufsize[0][Ey_->isDual(0)]-1;
        unsigned int iBz=istart[0][Bz_m->isDual(0)] + bufsize[0][Bz_m->isDual(0)]-1;
        unsigned int iEz=istart[0][Ez_->isDual(0)]  + bufsize[0][Ez_->isDual(0)]-1;
        unsigned int iBy=istart[0][By_m->isDual(0)] + bufsize[0][By_m->isDual(0)]-1;
        
        poynting_inst[1][0]=0.5*timestep*((*Ey_)(iEy) * ((*Bz_m)(iBz-1) + (*Bz_m)(iBz)) - 
                                          (*Ez_)(iEz) * ((*By_m)(iBy-1) + (*By_m)(iBy)));
        poynting[1][0] -= poynting_inst[1][0];
        
    }    
}

void ElectroMagn1D::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    Field1D* field1D=static_cast<Field1D*>(my_field);
    
    vector<double> pos(1);
    pos[0] = dx * ((double)(patch->getCellStartingGlobalIndex(0))+(field1D->isDual(0)?-0.5:0.));
    int N = (int)field1D->dims()[0];
    
    // USING UNSIGNED INT CREATES PB WITH PERIODIC BCs
    for (int i=0 ; i<N ; i++) {
        (*field1D)(i) += profile->valueAt(pos);
        pos[0] += dx;
    }
        
    if(emBoundCond[0]) emBoundCond[0]->save_fields_BC1D(my_field);
    if(emBoundCond[1]) emBoundCond[1]->save_fields_BC1D(my_field);
}


void ElectroMagn1D::initAntennas(Patch* patch)
{
    
    // Filling the space profiles of antennas
    for (unsigned int i=0; i<antennas.size(); i++) {
        if      (antennas[i].fieldName == "Jx")
            antennas[i].field = new Field1D(dimPrim, 0, false, "Jx");
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = new Field1D(dimPrim, 1, false, "Jy");
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = new Field1D(dimPrim, 2, false, "Jz");
        
        if (antennas[i].field) 
            applyExternalField(antennas[i].field, antennas[i].space_profile, patch);
    }

}
