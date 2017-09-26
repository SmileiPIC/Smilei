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
    for (unsigned int imode=0; imode<nmodes; imode++) {
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            ostringstream species_mode_name("");
            species_mode_name << vecSpecies[ispec]->species_type << "_mode_" << imode;
            Jx_s[imode*n_species+ispec]  = new cField2D(("Jx_" + species_mode_name.str()).c_str(), dimPrim);
            Jr_s[imode*n_species+ispec]  = new cField2D(("Jr_" + species_mode_name.str()).c_str(), dimPrim);
            Jt_s[imode*n_species+ispec]  = new cField2D(("Jt_" + species_mode_name.str()).c_str(), dimPrim);
            rho_s[imode*n_species+ispec] = new cField2D(("Rho_"+ species_mode_name.str()).c_str(), dimPrim);
        }
    }
    
}//END constructor Electromagn3D


ElectroMagn3DRZ::ElectroMagn3DRZ( ElectroMagn3DRZ* emFields, Params &params, Patch* patch ) : 
    ElectroMagn(emFields, params, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax())
{
    
    initElectroMagn3DRZQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int imode=0; imode<nmodes; imode++) {
        for (unsigned int ispec=0; ispec<n_species; ispec++) {

            int ifield = imode*n_species+ispec;

            if ( emFields->Jx_s[ifield] != NULL ) {
                if ( emFields->Jx_s[ifield]->data_ != NULL )
                    Jx_s[ifield]  = new cField2D(dimPrim, 0, false, emFields->Jx_s[ifield]->name);
                else
                    Jx_s[ifield]  = new cField2D(emFields->Jx_s[ifield]->name, dimPrim);
            }
            if ( emFields->Jr_s[ifield] != NULL ) {
                if ( emFields->Jr_s[ifield]->data_ != NULL )
                    Jr_s[ifield]  = new cField2D(dimPrim, 1, false, emFields->Jr_s[ifield]->name);
                else
                    Jr_s[ifield]  = new cField2D(emFields->Jr_s[ifield]->name, dimPrim);
            }
            if ( emFields->Jt_s[ifield] != NULL ) {
                if ( emFields->Jt_s[ifield]->data_ != NULL )
                    Jt_s[ifield]  = new cField2D(dimPrim, 2, false, emFields->Jt_s[ifield]->name);
                else
                    Jt_s[ifield]  = new cField2D(emFields->Jt_s[ifield]->name, dimPrim);
            }
            if ( emFields->rho_s[ifield] != NULL ) {
                if ( emFields->rho_s[ifield]->data_ != NULL )
                    rho_s[ifield] = new cField2D(dimPrim, emFields->rho_s[ifield]->name );
                else
                    rho_s[ifield]  = new cField2D(emFields->rho_s[ifield]->name, dimPrim);
            }
        }

    }


}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn3D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::initElectroMagn3DRZQuantities(Params &params, Patch* patch)
{

    // Species charge currents and density
    Jx_s.resize(n_species*nmodes);
    Jr_s.resize(n_species*nmodes);
    Jt_s.resize(n_species*nmodes);
    rho_s.resize(n_species*nmodes);
    for (unsigned int ispec=0; ispec<n_species*nmodes; ispec++) {
        Jx_s[ispec]  = NULL;
        Jr_s[ispec]  = NULL;
        Jt_s[ispec]  = NULL;
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

    Ex_.resize(nmodes);
    Er_.resize(nmodes);
    Et_.resize(nmodes);
    Bx_.resize(nmodes);
    Br_.resize(nmodes);
    Bt_.resize(nmodes);
    Bx_m.resize(nmodes);
    Br_m.resize(nmodes);
    Bt_m.resize(nmodes);
    
    // Total charge currents and densities
    Jx_.resize(nmodes);
    Jr_.resize(nmodes);
    Jt_.resize(nmodes);
    rho_.resize(nmodes);
    
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        ostringstream mode_id("");
        mode_id << "_mode_" << imode;

        Ex_[imode]  = new cField2D(dimPrim, 0, false, ("Ex"+mode_id.str()).c_str() );
        Er_[imode]  = new cField2D(dimPrim, 1, false, ("Er"+mode_id.str()).c_str() );
        Et_[imode]  = new cField2D(dimPrim, 2, false, ("Et"+mode_id.str()).c_str() );
        Bx_[imode]  = new cField2D(dimPrim, 0, true,  ("Bx"+mode_id.str()).c_str() );
        Br_[imode]  = new cField2D(dimPrim, 1, true,  ("Br"+mode_id.str()).c_str() );
        Bt_[imode]  = new cField2D(dimPrim, 2, true,  ("Bt"+mode_id.str()).c_str() );
        Bx_m[imode] = new cField2D(dimPrim, 0, true,  ("Bx_m"+mode_id.str()).c_str() );
        Br_m[imode] = new cField2D(dimPrim, 1, true,  ("Br_m"+mode_id.str()).c_str() );
        Bt_m[imode] = new cField2D(dimPrim, 2, true,  ("Bt_m"+mode_id.str()).c_str() );
    
        // Total charge currents and densities
        Jx_[imode]   = new cField2D(dimPrim, 0, false, ("Jx"+mode_id.str()).c_str() );
        Jr_[imode]   = new cField2D(dimPrim, 1, false, ("Jr"+mode_id.str()).c_str() );
        Jt_[imode]   = new cField2D(dimPrim, 2, false, ("Jt"+mode_id.str()).c_str() );
        rho_[imode]  = new cField2D(dimPrim, ("Rho"+mode_id.str()).c_str() );
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


void ElectroMagn3DRZ::finishInitialization(int nspecies, Patch* patch)
{
    // Fill allfields
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        allFields.push_back( Ex_[imode] );
        allFields.push_back( Er_[imode] );
        allFields.push_back( Et_[imode] );
        allFields.push_back( Bx_[imode] );
        allFields.push_back( Br_[imode] );
        allFields.push_back( Bt_[imode] );
        allFields.push_back( Bx_m[imode] );
        allFields.push_back( Br_m[imode] );
        allFields.push_back( Bt_m[imode] );
        allFields.push_back( Jx_[imode] );
        allFields.push_back( Jr_[imode] );
        allFields.push_back( Jt_[imode] );
        allFields.push_back( rho_[imode] );
    }

    for (unsigned int ispec=0; ispec<nspecies*nmodes; ispec++) {
        allFields.push_back(Jx_s[ispec] );
        allFields.push_back(Jr_s[ispec] );
        allFields.push_back(Jt_s[ispec] );
        allFields.push_back(rho_s[ispec]);
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn3DRZ
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3DRZ::~ElectroMagn3DRZ()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        delete Ex_[imode];
        delete Er_[imode];
        delete Et_[imode];
        delete Bx_[imode];
        delete Br_[imode];
        delete Bt_[imode];
        delete Bx_m[imode];
        delete Br_m[imode];
        delete Bt_m[imode];

        delete Jx_[imode];
        delete Jr_[imode];
        delete Jt_[imode];
        delete rho_[imode];
    }

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
    #ifdef _TODO_RZ
    cField2D* rho2D = static_cast<cField2D*>(rho_);

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

    phi_ = new cField2D(dimPrim);    // scalar potential
    r_   = new cField2D(dimPrim);    // residual vector
    p_   = new cField2D(dimPrim);    // direction vector
    Ap_  = new cField2D(dimPrim);    // A*p vector

    
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*phi_)(i,j)   = 0.0;
            (*r_)(i,j)     = -(*rho2D)(i,j);
            (*p_)(i,j)     = (*r_)(i,j);
        }//j
    }//i
    #endif

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
    #ifdef _TODO_RZ
    #endif
} // compute_pAp

double ElectroMagn3DRZ::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    #ifdef _TODO_RZ
    #endif
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
    #ifdef _TODO_RZ
    #endif

    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;

} // initE


void ElectroMagn3DRZ::centeringE( std::vector<double> E_Add )
{
    cField2D* Ex2D  = static_cast<cField2D*>(Ex_);
    cField2D* Er2D  = static_cast<cField2D*>(Er_);
    cField2D* Et2D  = static_cast<cField2D*>(Et_);

    // Centering electrostatic fields
    for (unsigned int i=0; i<nx_d; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Ex2D)(i,j) += E_Add[0];
        }
    }
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_d; j++) {
            (*Er2D)(i,j) += E_Add[1];
        }
    }
    for (unsigned int i=0; i<nx_p; i++) {
        for (unsigned int j=0; j<ny_p; j++) {
            (*Et2D)(i,j) += E_Add[2];
        }
    }
    #ifdef _TODO_RZ
    #endif

} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods 
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::saveMagneticFields()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        // Static cast of the fields
        cField2D* Bx3DRZ    = static_cast<cField2D*>(Bx_[imode]);
        cField2D* Br3DRZ    = static_cast<cField2D*>(Br_[imode]);
        cField2D* Bt3DRZ    = static_cast<cField2D*>(Bt_[imode]);
        cField2D* Bx3DRZ_m = static_cast<cField2D*>(Bx_m[imode]);
        cField2D* Br3DRZ_m = static_cast<cField2D*>(Br_m[imode]);
        cField2D* Bt3DRZ_m = static_cast<cField2D*>(Bt_m[imode]);
    
        // Magnetic field Bx^(p,d)
        memcpy(&((*Bx3DRZ_m)(0,0)), &((*Bx3DRZ)(0,0)),nx_p*ny_d*sizeof(complex<double>) );
    
        // Magnetic field Br^(d,p)
        memcpy(&((*Br3DRZ_m)(0,0)), &((*Br3DRZ)(0,0)),nx_d*ny_p*sizeof(complex<double>) );
    
        // Magnetic field Bt^(d,d)
        memcpy(&((*Bt3DRZ_m)(0,0)), &((*Bt3DRZ)(0,0)),nx_d*ny_d*sizeof(complex<double>) );
    }

}//END saveMagneticFields


// Create a new field
Field * ElectroMagn3DRZ::createField(string fieldname)
{
    if     (fieldname.substr(0,2)=="Ex" ) return new cField2D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Er" ) return new cField2D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Et" ) return new cField2D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,2)=="Bx" ) return new cField2D(dimPrim, 0, true,  fieldname);
    else if(fieldname.substr(0,2)=="Br" ) return new cField2D(dimPrim, 1, true,  fieldname);
    else if(fieldname.substr(0,2)=="Bt" ) return new cField2D(dimPrim, 2, true,  fieldname);
    else if(fieldname.substr(0,2)=="Jx" ) return new cField2D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,2)=="Jr" ) return new cField2D(dimPrim, 1, false, fieldname);
    else if(fieldname.substr(0,2)=="Jt" ) return new cField2D(dimPrim, 2, false, fieldname);
    else if(fieldname.substr(0,3)=="Rho") return new cField2D(dimPrim, fieldname );
    
    ERROR("Cannot create field "<<fieldname);
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::centerMagneticFields()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {

        // Static cast of the fields
        cField2D* Bx3DRZ    = static_cast<cField2D*>(Bx_[imode]);
        cField2D* Br3DRZ    = static_cast<cField2D*>(Br_[imode]);
        cField2D* Bt3DRZ    = static_cast<cField2D*>(Bt_[imode]);
        cField2D* Bx3DRZ_m  = static_cast<cField2D*>(Bx_m[imode]);
        cField2D* Br3DRZ_m  = static_cast<cField2D*>(Br_m[imode]);
        cField2D* Bt3DRZ_m  = static_cast<cField2D*>(Bt_m[imode]);
    
        // Magnetic field Bx^(p,d)
        for (unsigned int i=0 ; i<nx_p ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bx3DRZ_m)(i,j) = ( (*Bx3DRZ)(i,j) + (*Bx3DRZ_m)(i,j) )*0.5;
            }
        }
    
        // Magnetic field Br^(d,p)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_p ; j++) {
                (*Br3DRZ_m)(i,j) = ( (*Br3DRZ)(i,j) + (*Br3DRZ_m)(i,j) )*0.5;
            }
        }
    
        // Magnetic field Bt^(d,d)
        for (unsigned int i=0 ; i<nx_d ; i++) {
            for (unsigned int j=0 ; j<ny_d ; j++) {
                (*Bt3DRZ_m)(i,j) = ( (*Bt3DRZ)(i,j) + (*Bt3DRZ_m)(i,j) )*0.5;
            } // end for j
        } // end for i

    }
    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::binomialCurrentFilter()
{
    ERROR("Binomial current filtering not yet implemented in 3DRZ");
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computeTotalRhoJ()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {

        // static cast of the total currents and densities
        cField2D* Jx3DRZ    = static_cast<cField2D*>(Jx_[imode]);
        cField2D* Jr3DRZ    = static_cast<cField2D*>(Jr_[imode]);
        cField2D* Jt3DRZ    = static_cast<cField2D*>(Jt_[imode]);
        cField2D* rho3DRZ   = static_cast<cField2D*>(rho_[imode]);    
    
        // -----------------------------------
        // Species currents and charge density
        // -----------------------------------
        for (unsigned int ispec=0; ispec<n_species; ispec++) {

            int ifield = imode*n_species+ispec;

            if( Jx_s[ifield] ) {
                cField2D* Jx2D_s  = static_cast<cField2D*>(Jx_s[ifield]);
                for (unsigned int i=0 ; i<=nx_p ; i++)
                    for (unsigned int j=0 ; j<ny_p ; j++)
                        (*Jx3DRZ)(i,j) += (*Jx2D_s)(i,j);
            }
            if( Jr_s[ifield] ) {
                cField2D* Jr2D_s  = static_cast<cField2D*>(Jr_s[ifield]);
                for (unsigned int i=0 ; i<nx_p ; i++)
                    for (unsigned int j=0 ; j<=ny_p ; j++)
                        (*Jr3DRZ)(i,j) += (*Jr2D_s)(i,j);
            }
            if( Jt_s[ifield] ) {
                cField2D* Jt2D_s  = static_cast<cField2D*>(Jt_s[ifield]);
                for (unsigned int i=0 ; i<nx_p ; i++)
                    for (unsigned int j=0 ; j<ny_p ; j++)
                        (*Jt3DRZ)(i,j) += (*Jt2D_s)(i,j);
            }
            if( rho_s[ifield] ) {
                cField2D* rho2D_s  = static_cast<cField2D*>(rho_s[ifield]);
                for (unsigned int i=0 ; i<nx_p ; i++)
                    for (unsigned int j=0 ; j<ny_p ; j++)
                        (*rho3DRZ)(i,j) += (*rho2D_s)(i,j);
            }
        
        }//END loop on species ispec
        
    }//END loop on mmodes
    
} //END computeTotalRhoJ


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computePoynting() {

    cField2D* Ex2D     = static_cast<cField2D*>(Ex_);
    cField2D* Er2D     = static_cast<cField2D*>(Er_);
    cField2D* Et2D     = static_cast<cField2D*>(Et_);
    cField2D* Bx2D_m   = static_cast<cField2D*>(Bx_m);
    cField2D* Br2D_m   = static_cast<cField2D*>(Br_m);
    cField2D* Bt2D_m   = static_cast<cField2D*>(Bt_m);

    if (isXmin) {
        unsigned int iEr=istart[0][Er2D->isDual(0)];
        unsigned int iBt=istart[0][Bt2D_m->isDual(0)];
        unsigned int iEt=istart[0][Et2D->isDual(0)];
        unsigned int iBr=istart[0][Br2D_m->isDual(0)];
        
        unsigned int jEr=istart[1][Er2D->isDual(1)];
        unsigned int jBt=istart[1][Bt2D_m->isDual(1)];
        unsigned int jEt=istart[1][Et2D->isDual(1)];
        unsigned int jBr=istart[1][Br2D_m->isDual(1)];
        

        for (unsigned int j=0; j<=bufsize[1][Et2D->isDual(1)]; j++) {
            #ifdef _TODO_RZ            
            double Er__ = 0.5*((*Er2D)(iEr,jEr+j) + (*Er2D)(iEr, jEr+j+1));
            double Bt__ = 0.25*((*Bt2D_m)(iBt,jBt+j)+(*Bt2D_m)(iBt+1,jBt+j)+(*Bt2D_m)(iBt,jBt+j+1)+(*Bt2D_m)(iBt+1,jBt+j+1));
            double Et__ = (*Et2D)(iEt,jEt+j);
            double Br__ = 0.5*((*Br2D_m)(iBr,jBr+j) + (*Br2D_m)(iBr+1, jBr+j));

            poynting_inst[0][0] = dy*timestep*(Er__*Bt__ - Et__*Br__);
            #endif
            poynting[0][0]+= poynting_inst[0][0];

        }
        
    }//if Xmin
    
    
    if (isXmax) {
        
        unsigned int iEr=istart[0][Er2D->isDual(0)]  + bufsize[0][Er2D->isDual(0)] -1;
        unsigned int iBt=istart[0][Bt2D_m->isDual(0)] + bufsize[0][Bt2D_m->isDual(0)]-1;
        unsigned int iEt=istart[0][Et2D->isDual(0)]  + bufsize[0][Et2D->isDual(0)] -1;
        unsigned int iBr=istart[0][Br2D_m->isDual(0)] + bufsize[0][Br2D_m->isDual(0)]-1;
        
        unsigned int jEr=istart[1][Er2D->isDual(1)];
        unsigned int jBt=istart[1][Bt2D_m->isDual(1)];
        unsigned int jEt=istart[1][Et2D->isDual(1)];
        unsigned int jBr=istart[1][Br2D_m->isDual(1)];
        
        for (unsigned int j=0; j<=bufsize[1][Et2D->isDual(1)]; j++) {
            #ifdef _TODO_RZ            
          
            double Er__ = 0.5*((*Er2D)(iEr,jEr+j) + (*Er2D)(iEr, jEr+j+1));
            double Bt__ = 0.25*((*Bt2D_m)(iBt,jBt+j)+(*Bt2D_m)(iBt+1,jBt+j)+(*Bt2D_m)(iBt,jBt+j+1)+(*Bt2D_m)(iBt+1,jBt+j+1));
            double Et__ = (*Et2D)(iEt,jEt+j);
            double Br__ = 0.5*((*Br2D_m)(iBr,jBr+j) + (*Br2D_m)(iBr+1, jBr+j));
            
            poynting_inst[1][0] = dy*timestep*(Er__*Bt__ - Et__*Br__);
            #endif
            poynting[1][0]+= poynting_inst[1][0];

        }
        
    }//if Xmax
    
    if (isYmin) {
        
        unsigned int iEt=istart[0][Et_->isDual(0)];
        unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
        unsigned int iEx=istart[0][Ex_->isDual(0)];
        unsigned int iBt=istart[0][Bt_m->isDual(0)]; 
        
        unsigned int jEt=istart[1][Et_->isDual(1)];
        unsigned int jBx=istart[1][Bx_m->isDual(1)];
        unsigned int jEx=istart[1][Ex_->isDual(1)];
        unsigned int jBt=istart[1][Bt_m->isDual(1)];

        for (unsigned int i=0; i<=bufsize[0][Et2D->isDual(0)]; i++) {
            #ifdef _TODO_RZ            
            double Et__ = (*Et2D)(iEt+i,jEt);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bt__ = 0.25*((*Bt2D_m)(iBt+i,jBt)+(*Bt2D_m)(iBt+i+1,jBt)+(*Bt2D_m)(iBt+i,jBt+1)+(*Bt2D_m)(iBt+i+1,jBt+1));
            
            poynting_inst[0][1] = dx*timestep*(Et__*Bx__ - Ex__*Bt__);
            #endif
            poynting[0][1] += poynting_inst[0][1];
        }

    }// if Ymin
    
    if (isYmax) {

        unsigned int iEt=istart[0][Et2D->isDual(0)];
        unsigned int iBx=istart[0][Bx2D_m->isDual(0)];
        unsigned int iEx=istart[0][Ex2D->isDual(0)];
        unsigned int iBt=istart[0][Bt2D_m->isDual(0)];
        
        unsigned int jEt=istart[1][Et2D->isDual(1)]  + bufsize[1][Et2D->isDual(1)] -1;
        unsigned int jBx=istart[1][Bx2D_m->isDual(1)] + bufsize[1][Bx2D_m->isDual(1)]-1;
        unsigned int jEx=istart[1][Ex2D->isDual(1)]  + bufsize[1][Ex2D->isDual(1)] -1;
        unsigned int jBt=istart[1][Bt2D_m->isDual(1)] + bufsize[1][Bt2D_m->isDual(1)]-1;
        
        for (unsigned int i=0; i<=bufsize[0][Et2D->isDual(0)]; i++) {
            #ifdef _TODO_RZ            
            double Et__ = (*Et2D)(iEt+i,jEt);
            double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
            double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
            double Bt__ = 0.25*((*Bt2D_m)(iBt+i,jBt)+(*Bt2D_m)(iBt+i+1,jBt)+(*Bt2D_m)(iBt+i,jBt+1)+(*Bt2D_m)(iBt+i+1,jBt+1));
            
            poynting_inst[1][1] = dx*timestep*(Et__*Bx__ - Ex__*Bt__);
            #endif
            poynting[1][1] += poynting_inst[1][1];
        }

    }//if Ymax

}

void ElectroMagn3DRZ::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    
    cField2D* field2D=static_cast<cField2D*>(my_field);
    
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
            antennas[i].field = new cField2D(dimPrim, 0, false, "Jx");
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = new cField2D(dimPrim, 1, false, "Jy");
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = new cField2D(dimPrim, 2, false, "Jz");
        else {
            ERROR("Antenna cannot be applied to field "<<antennas[i].fieldName);
        }
        
        if (antennas[i].field) 
            applyExternalField(antennas[i].field, antennas[i].space_profile, patch);
    }

}

