#include "ElectroMagn3DRZ.h"

#include <cmath>

#include <iostream>
#include <sstream>
#include <complex>
#include "dcomplex.h"
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
ElectroMagn3DRZ::ElectroMagn3DRZ(Params &params, DomainDecomposition* domain_decomposition, vector<Species*>& vecSpecies, Patch* patch) : 
  ElectroMagn(params, domain_decomposition, vecSpecies, patch),
isYmin(patch->isYmin()),
isYmax(patch->isYmax())
{    
    
    initElectroMagn3DRZQuantities(params, patch);
    
    // Charge currents currents and density for each species
    for (unsigned int imode=0; imode<nmodes; imode++) {
        for (unsigned int ispec=0; ispec<n_species; ispec++) {
            ostringstream species_mode_name("");
            species_mode_name << vecSpecies[ispec]->name << "_mode_" << imode;
            Jl_s[imode*n_species+ispec]  = new cField2D(("Jx_" + species_mode_name.str()).c_str(), dimPrim);
            Jr_s[imode*n_species+ispec]  = new cField2D(("Jr_" + species_mode_name.str()).c_str(), dimPrim);
            Jt_s[imode*n_species+ispec]  = new cField2D(("Jt_" + species_mode_name.str()).c_str(), dimPrim);
            rho_RZ_s[imode*n_species+ispec] = new cField2D(("Rho_"+ species_mode_name.str()).c_str(), dimPrim);
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

            if ( emFields->Jl_s[ifield] != NULL ) {
                if ( emFields->Jl_s[ifield]->cdata_ != NULL )
                    Jl_s[ifield]  = new cField2D(dimPrim, 0, false, emFields->Jl_s[ifield]->name);
                else
                    Jl_s[ifield]  = new cField2D(emFields->Jl_s[ifield]->name, dimPrim);
            }
            if ( emFields->Jr_s[ifield] != NULL ) {
                if ( emFields->Jr_s[ifield]->cdata_ != NULL )
                    Jr_s[ifield]  = new cField2D(dimPrim, 1, false, emFields->Jr_s[ifield]->name);
                else
                    Jr_s[ifield]  = new cField2D(emFields->Jr_s[ifield]->name, dimPrim);
            }
            if ( emFields->Jt_s[ifield] != NULL ) {
                if ( emFields->Jt_s[ifield]->cdata_ != NULL )
                    Jt_s[ifield]  = new cField2D(dimPrim, 2, false, emFields->Jt_s[ifield]->name);
                else
                    Jt_s[ifield]  = new cField2D(emFields->Jt_s[ifield]->name, dimPrim);
            }
            if ( emFields->rho_RZ_s[ifield] != NULL ) {
                if ( emFields->rho_RZ_s[ifield]->cdata_ != NULL )
                    rho_RZ_s[ifield] = new cField2D(dimPrim, emFields->rho_RZ_s[ifield]->name );
                else
                    rho_RZ_s[ifield]  = new cField2D(emFields->rho_RZ_s[ifield]->name, dimPrim);
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
    Jl_s.resize(n_species*nmodes);
    Jr_s.resize(n_species*nmodes);
    Jt_s.resize(n_species*nmodes);
    rho_RZ_s.resize(n_species*nmodes);
    for (unsigned int ispec=0; ispec<n_species*nmodes; ispec++) {
        Jl_s[ispec]  = NULL;
        Jr_s[ispec]  = NULL;
        Jt_s[ispec]  = NULL;
        rho_RZ_s[ispec] = NULL;
    }

    // --------------------------------------------------
    // Calculate quantities related to the simulation box
    // --------------------------------------------------
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dl       = cell_length[0];
    dt_ov_dl = timestep/dl;
    dl_ov_dt = 1.0/dt_ov_dl;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dr       = cell_length[1];
    dt_ov_dr = timestep/dr;
    dr_ov_dt = 1.0/dt_ov_dr;
    j_glob_ = patch->getCellStartingGlobalIndex(1);

    
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
    nl_p = n_space[0]+1+2*oversize[0];
    nl_d = n_space[0]+2+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    nr_p = n_space[1]+1+2*oversize[1];
    nr_d = n_space[1]+2+2*oversize[1];
    
    // Allocation of the EM fields

    El_.resize(nmodes);
    Er_.resize(nmodes);
    Et_.resize(nmodes);
    Bl_.resize(nmodes);
    Br_.resize(nmodes);
    Bt_.resize(nmodes);
    Bl_m.resize(nmodes);
    Br_m.resize(nmodes);
    Bt_m.resize(nmodes);
    
    // Total charge currents and densities
    Jl_.resize(nmodes);
    Jr_.resize(nmodes);
    Jt_.resize(nmodes);
    rho_RZ_.resize(nmodes);
    
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        ostringstream mode_id("");
        mode_id << "_mode_" << imode;

        El_[imode]  = new cField2D(dimPrim, 0, false, ("Ex"+mode_id.str()).c_str() );
        Er_[imode]  = new cField2D(dimPrim, 1, false, ("Er"+mode_id.str()).c_str() );
        Et_[imode]  = new cField2D(dimPrim, 2, false, ("Et"+mode_id.str()).c_str() );
        Bl_[imode]  = new cField2D(dimPrim, 0, true,  ("Bx"+mode_id.str()).c_str() );
        Br_[imode]  = new cField2D(dimPrim, 1, true,  ("Br"+mode_id.str()).c_str() );
        Bt_[imode]  = new cField2D(dimPrim, 2, true,  ("Bt"+mode_id.str()).c_str() );
        Bl_m[imode] = new cField2D(dimPrim, 0, true,  ("Bx_m"+mode_id.str()).c_str() );
        Br_m[imode] = new cField2D(dimPrim, 1, true,  ("Br_m"+mode_id.str()).c_str() );
        Bt_m[imode] = new cField2D(dimPrim, 2, true,  ("Bt_m"+mode_id.str()).c_str() );
    
        // Total charge currents and densities
        Jl_[imode]   = new cField2D(dimPrim, 0, false, ("Jx"+mode_id.str()).c_str() );
        Jr_[imode]   = new cField2D(dimPrim, 1, false, ("Jr"+mode_id.str()).c_str() );
        Jt_[imode]   = new cField2D(dimPrim, 2, false, ("Jt"+mode_id.str()).c_str() );
        rho_RZ_[imode]  = new cField2D(dimPrim, ("Rho"+mode_id.str()).c_str() );
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
     MESSAGE("index_bc_min / index_bc_max / nl_p / nl_d" << index_bc_min[0]
     << " " << index_bc_max[0] << " " << nl_p<< " " << nl_d);
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
        allFields.push_back( El_[imode] );
        allFields.push_back( Er_[imode] );
        allFields.push_back( Et_[imode] );
        allFields.push_back( Bl_[imode] );
        allFields.push_back( Br_[imode] );
        allFields.push_back( Bt_[imode] );
        allFields.push_back( Bl_m[imode] );
        allFields.push_back( Br_m[imode] );
        allFields.push_back( Bt_m[imode] );
        allFields.push_back( Jl_[imode] );
        allFields.push_back( Jr_[imode] );
        allFields.push_back( Jt_[imode] );
        allFields.push_back( rho_RZ_[imode] );
    }

    for (int ispec=0; ispec<nspecies*(int)nmodes; ispec++) {
        allFields.push_back(Jl_s[ispec] );
        allFields.push_back(Jr_s[ispec] );
        allFields.push_back(Jt_s[ispec] );
        allFields.push_back(rho_RZ_s[ispec]);
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Electromagn3DRZ
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3DRZ::~ElectroMagn3DRZ()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        delete El_[imode];
        delete Er_[imode];
        delete Et_[imode];
        delete Bl_[imode];
        delete Br_[imode];
        delete Bt_[imode];
        delete Bl_m[imode];
        delete Br_m[imode];
        delete Bt_m[imode];

        delete Jl_[imode];
        delete Jr_[imode];
        delete Jt_[imode];
        delete rho_RZ_[imode];
    }

}//END ElectroMagn3DRZ


void ElectroMagn3DRZ::restartRhoJ()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        Jl_[imode] ->put_to(0.);
        Jr_[imode] ->put_to(0.);
        Jt_[imode] ->put_to(0.);
        rho_RZ_[imode]->put_to(0.);
    }
}

void ElectroMagn3DRZ::restartRhoJs()
{
    for (unsigned int ispec=0 ; ispec < n_species*nmodes ; ispec++) {
        if( Jl_s [ispec] ) Jl_s [ispec]->put_to(0.);
        if( Jr_s [ispec] ) Jr_s [ispec]->put_to(0.);
        if( Jt_s [ispec] ) Jt_s [ispec]->put_to(0.);
        if( rho_RZ_s[ispec] ) rho_RZ_s[ispec]->put_to(0.);
    }
    
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        Jl_[imode] ->put_to(0.);
        Jr_ [imode]->put_to(0.);
        Jt_ [imode]->put_to(0.);
        rho_RZ_[imode]->put_to(0.);
    }
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
    index_max_p_[0] = nl_p - 2 - oversize[0];
    index_max_p_[1] = nr_p - 2 - oversize[1];
    if (patch->isXmin()) {
        index_min_p_[0] = 0;
    }
    if (patch->isXmax()) {
        index_max_p_[0] = nl_p-1;
    }

    phi_ = new cField2D(dimPrim);    // scalar potential
    r_   = new cField2D(dimPrim);    // residual vector
    p_   = new cField2D(dimPrim);    // direction vector
    Ap_  = new cField2D(dimPrim);    // A*p vector

    
    for (unsigned int i=0; i<nl_p; i++) {
        for (unsigned int j=0; j<nr_p; j++) {
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
    for(unsigned int i=0; i<nl_p; i++) {
        for(unsigned int j=0; j<nr_p; j++) {
            (*phi_)(i,j) += alpha_k * (*p_)(i,j);
            (*r_)(i,j)   -= alpha_k * (*Ap_)(i,j);
        }
    }

} // update_pand_r

void ElectroMagn3DRZ::update_p(double rnew_dot_rnew, double r_dot_r)
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for (unsigned int i=0; i<nl_p; i++) {
        for(unsigned int j=0; j<nr_p; j++) {
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
    cField2D* Ey2D  = static_cast<cField2D*>(Ey_);
    cField2D* Ez2D  = static_cast<cField2D*>(Ez_);

    // Centering electrostatic fields
    for (unsigned int i=0; i<nl_d; i++) {
        for (unsigned int j=0; j<nr_p; j++) {
            (*Ex2D)(i,j) += E_Add[0];
        }
    }
    for (unsigned int i=0; i<nl_p; i++) {
        for (unsigned int j=0; j<nr_d; j++) {
            (*Ey2D)(i,j) += E_Add[1];
        }
    }
    for (unsigned int i=0; i<nl_p; i++) {
        for (unsigned int j=0; j<nr_p; j++) {
            (*Ez2D)(i,j) += E_Add[2];
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
void ElectroMagn3DRZ::saveMagneticFields(bool is_spectral)
{
    if (is_spectral)
        ERROR("Not implemented");
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        // Static cast of the fields
        cField2D* Bl3DRZ    = static_cast<cField2D*>(Bl_[imode]);
        cField2D* Br3DRZ    = static_cast<cField2D*>(Br_[imode]);
        cField2D* Bt3DRZ    = static_cast<cField2D*>(Bt_[imode]);
        cField2D* Bl3D_RZ_m = static_cast<cField2D*>(Bl_m[imode]);
        cField2D* Br3D_RZ_m = static_cast<cField2D*>(Br_m[imode]);
        cField2D* Bt3D_RZ_m = static_cast<cField2D*>(Bt_m[imode]);
    
        // Magnetic field Bx^(p,d)
        memcpy(&((*Bl3D_RZ_m)(0,0)), &((*Bl3DRZ)(0,0)),nl_p*nr_d*sizeof(complex<double>) );
    
        // Magnetic field Br^(d,p)
        memcpy(&((*Br3D_RZ_m)(0,0)), &((*Br3DRZ)(0,0)),nl_d*nr_p*sizeof(complex<double>) );
    
        // Magnetic field Bt^(d,d)
        memcpy(&((*Bt3D_RZ_m)(0,0)), &((*Bt3DRZ)(0,0)),nl_d*nr_d*sizeof(complex<double>) );
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
    return NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::centerMagneticFields()
{
    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {

        // Static cast of the fields
        cField2D* Bl3DRZ    = static_cast<cField2D*>(Bl_[imode]);
        cField2D* Br3DRZ    = static_cast<cField2D*>(Br_[imode]);
        cField2D* Bt3DRZ    = static_cast<cField2D*>(Bt_[imode]);
        cField2D* Bl3D_RZ_m = static_cast<cField2D*>(Bl_m[imode]);
        cField2D* Br3D_RZ_m = static_cast<cField2D*>(Br_m[imode]);
        cField2D* Bt3D_RZ_m = static_cast<cField2D*>(Bt_m[imode]);
    
        // Magnetic field Bx^(p,d,d)
        for (unsigned int i=0 ; i<nl_p ; i++) {
            for (unsigned int j=0 ; j<nr_d ; j++) {
                (*Bl3D_RZ_m)(i,j) = ( (*Bl3DRZ)(i,j) + (*Bl3D_RZ_m)(i,j) )*0.5;
            }
        }
    
        // Magnetic field Br^(d,p,d)
        for (unsigned int i=0 ; i<nl_d ; i++) {
            for (unsigned int j=0 ; j<nr_p ; j++) {
                (*Br3D_RZ_m)(i,j) = ( (*Br3DRZ)(i,j) + (*Br3D_RZ_m)(i,j) )*0.5;
            }
        }
    
        // Magnetic field Bt^(d,d,p)
        for (unsigned int i=0 ; i<nl_d ; i++) {
            for (unsigned int j=0 ; j<nr_d ; j++) {
                (*Bt3D_RZ_m)(i,j) = ( (*Bt3DRZ)(i,j) + (*Bt3D_RZ_m)(i,j) )*0.5;
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
        cField2D* JlRZ    = static_cast<cField2D*>(Jl_[imode]);
        cField2D* JrRZ    = static_cast<cField2D*>(Jr_[imode]);
        cField2D* JtRZ    = static_cast<cField2D*>(Jt_[imode]);
        cField2D* rhoRZ   = static_cast<cField2D*>(rho_RZ_[imode]);    
        //MESSAGE("c");
        // -----------------------------------
        // Species currents and charge density
        // -----------------------------------
        for (unsigned int ispec=0; ispec<n_species; ispec++) {

            int ifield = imode*n_species+ispec;
           // MESSAGE("cc");
	   // MESSAGE(Jl_s.size());
	   // MESSAGE(Jl_.size());
	   // MESSAGE(ifield);
	    if( Jl_s[ifield] ) {
                cField2D* Jl2D_s  = static_cast<cField2D*>(Jl_s[ifield]);
                MESSAGE(Jl2D_s->dims_[0]);
	        MESSAGE(Jl2D_s->dims_[1]);
		MESSAGE(JlRZ->dims_[0]);
		MESSAGE(JlRZ->dims_[1]);
                for (unsigned int i=0 ; i<=nl_p ; i++){
		    //MESSAGE("here");
		    //MESSAGE(nr_p);
		    //MESSAGE(nl_p);
                    for (unsigned int j=0 ; j<nr_p ; j++){
			//MESSAGE("here i=" <<i << "  j="<<j);
                        (*JlRZ)(i,j) += (*Jl2D_s)(i,j);}}
            }
	    //MESSAGE("or here");
            if( Jr_s[ifield] ) {
                cField2D* Jr2D_s  = static_cast<cField2D*>(Jr_s[ifield]);
                for (unsigned int i=0 ; i<nl_p ; i++)
                    for (unsigned int j=0 ; j<=nr_p ; j++)
                        (*JrRZ)(i,j) += (*Jr2D_s)(i,j);
            }
            if( Jt_s[ifield] ) {
                cField2D* Jt2D_s  = static_cast<cField2D*>(Jt_s[ifield]);
                for (unsigned int i=0 ; i<nl_p ; i++)
                    for (unsigned int j=0 ; j<nr_p ; j++)
                        (*JtRZ)(i,j) += (*Jt2D_s)(i,j);
            }
            if( rho_RZ_s[ifield] ) {
                cField2D* rho2D_s  = static_cast<cField2D*>(rho_RZ_s[ifield]);
                for (unsigned int i=0 ; i<nl_p ; i++)
                    for (unsigned int j=0 ; j<nr_p ; j++)
                        (*rhoRZ)(i,j) += (*rho2D_s)(i,j);
            }
        
        }//END loop on species ispec
        
    }//END loop on mmodes
    //MESSAGE("totalRj");
} //END computeTotalRhoJ

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total susceptibility from species susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computeTotalEnvChi()
{ } //END computeTotalEnvChi


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3DRZ::computePoynting() {

    //cField2D* Ex2D     = static_cast<cField2D*>(Ex_);
    //cField2D* Ey2D     = static_cast<cField2D*>(Ey_);
    //cField2D* Ez2D     = static_cast<cField2D*>(Ez_);
    //cField2D* Bx2D_m   = static_cast<cField2D*>(Bx_m);
    //cField2D* By2D_m   = static_cast<cField2D*>(By_m);
    //cField2D* Bz2D_m   = static_cast<cField2D*>(Bz_m);

    //if (isXmin) {
    //    unsigned int iEy=istart[0][Ey2D->isDual(0)];
    //    unsigned int iBz=istart[0][Bz2D_m->isDual(0)];
    //    unsigned int iEz=istart[0][Ez2D->isDual(0)];
    //    unsigned int iBy=istart[0][By2D_m->isDual(0)];
    //    
    //    unsigned int jEy=istart[1][Ey2D->isDual(1)];
    //    unsigned int jBz=istart[1][Bz2D_m->isDual(1)];
    //    unsigned int jEz=istart[1][Ez2D->isDual(1)];
    //    unsigned int jBy=istart[1][By2D_m->isDual(1)];
    //    

    //    for (unsigned int j=0; j<=bufsize[1][Ez2D->isDual(1)]; j++) {
    //        #ifdef _TODO_RZ            
    //        double Ey__ = 0.5*((*Ey2D)(iEr,jEy+j) + (*Ey2D)(iEy, jEy+j+1));
    //        double Bz__ = 0.25*((*Bz2D_m)(iBz,jBz+j)+(*Bz2D_m)(iBz+1,jBz+j)+(*Bz2D_m)(iBz,jBz+j+1)+(*Bz2D_m)(iBz+1,jBz+j+1));
    //        double Ez__ = (*Ez2D)(iEz,jEz+j);
    //        double By__ = 0.5*((*By2D_m)(iBy,jBy+j) + (*By2D_m)(iBy+1, jBy+j));

    //        poynting_inst[0][0] = dr*timestep*(Ey__*Bz__ - Ez__*By__);
    //        #endif
    //        poynting[0][0]+= poynting_inst[0][0];

    //    }
    //    
    //}//if Xmin
    //
    //
    //if (isXmax) {
    //    
    //    unsigned int iEy=istart[0][Ey2D->isDual(0)]  + bufsize[0][Ey2D->isDual(0)] -1;
    //    unsigned int iBz=istart[0][Bz2D_m->isDual(0)] + bufsize[0][Bz2D_m->isDual(0)]-1;
    //    unsigned int iEz=istart[0][Ez2D->isDual(0)]  + bufsize[0][Ez2D->isDual(0)] -1;
    //    unsigned int iBy=istart[0][By2D_m->isDual(0)] + bufsize[0][By2D_m->isDual(0)]-1;
    //    
    //    unsigned int jEy=istart[1][Ey2D->isDual(1)];
    //    unsigned int jBz=istart[1][Bz2D_m->isDual(1)];
    //    unsigned int jEz=istart[1][Ez2D->isDual(1)];
    //    unsigned int jBy=istart[1][By2D_m->isDual(1)];
    //    
    //    for (unsigned int j=0; j<=bufsize[1][Ez2D->isDual(1)]; j++) {
    //        #ifdef _TODO_RZ            
    //      
    //        double Ey__ = 0.5*((*Ey2D)(iEy,jEy+j) + (*Ey2D)(iEr, jEy+j+1));
    //        double Bz__ = 0.25*((*Bz2D_m)(iBz,jBz+j)+(*Bz2D_m)(iBz+1,jBz+j)+(*Bz2D_m)(iBz,jBz+j+1)+(*Bz2D_m)(iBz+1,jBz+j+1));
    //        double Ez__ = (*Ez2D)(iEz,jEz+j);
    //        double By__ = 0.5*((*By2D_m)(iBy,jBy+j) + (*By2D_m)(iBy+1, jBy+j));
    //        
    //        poynting_inst[1][0] = dr*timestep*(Ey__*Bz__ - Ez__*By__);
    //        #endif
    //        poynting[1][0]+= poynting_inst[1][0];

    //    }
    //    
    //}//if Xmax
    //
    //if (isYmin) {
    //    
    //    unsigned int iEz=istart[0][Ez_->isDual(0)];
    //    unsigned int iBx=istart[0][Bx_m->isDual(0)]; 
    //    unsigned int iEx=istart[0][Ex_->isDual(0)];
    //    unsigned int iBz=istart[0][Bz_m->isDual(0)]; 
    //    
    //    unsigned int jEz=istart[1][Ez_->isDual(1)];
    //    unsigned int jBx=istart[1][Bx_m->isDual(1)];
    //    unsigned int jEx=istart[1][Ex_->isDual(1)];
    //    unsigned int jBz=istart[1][Bz_m->isDual(1)];

    //    for (unsigned int i=0; i<=bufsize[0][Ez2D->isDual(0)]; i++) {
    //        #ifdef _TODO_RZ            
    //        double Ez__ = (*Ez2D)(iEz+i,jEz);
    //        double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
    //        double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
    //        double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
    //        
    //        poynting_inst[0][1] = dl*timestep*(Ez__*Bx__ - Ex__*Bz__);
    //        #endif
    //        poynting[0][1] += poynting_inst[0][1];
    //    }

    //}// if Ymin
    //
    //if (isYmax) {

    //    unsigned int iEz=istart[0][Ez2D->isDual(0)];
    //    unsigned int iBx=istart[0][Bx2D_m->isDual(0)];
    //    unsigned int iEx=istart[0][Ex2D->isDual(0)];
    //    unsigned int iBz=istart[0][Bz2D_m->isDual(0)];
    //    
    //    unsigned int jEz=istart[1][Ez2D->isDual(1)]  + bufsize[1][Ez2D->isDual(1)] -1;
    //    unsigned int jBx=istart[1][Bx2D_m->isDual(1)] + bufsize[1][Bx2D_m->isDual(1)]-1;
    //    unsigned int jEx=istart[1][Ex2D->isDual(1)]  + bufsize[1][Ex2D->isDual(1)] -1;
    //    unsigned int jBz=istart[1][Bz2D_m->isDual(1)] + bufsize[1][Bz2D_m->isDual(1)]-1;
    //    
    //    for (unsigned int i=0; i<=bufsize[0][Ez2D->isDual(0)]; i++) {
    //        #ifdef _TODO_RZ            
    //        double Ez__ = (*Ez2D)(iEz+i,jEz);
    //        double Bx__ = 0.5*((*Bx2D_m)(iBx+i,jBx) + (*Bx2D_m)(iBx+i, jBx+1));
    //        double Ex__ = 0.5*((*Ex2D)(iEx+i,jEx) + (*Ex2D)(iEx+i+1, jEx));
    //        double Bz__ = 0.25*((*Bz2D_m)(iBz+i,jBz)+(*Bz2D_m)(iBz+i+1,jBz)+(*Bz2D_m)(iBz+i,jBz+1)+(*Bz2D_m)(iBz+i+1,jBz+1));
    //        
    //        poynting_inst[1][1] = dl*timestep*(Ez__*Bx__ - Ex__*Bz__);
    //        #endif
    //        poynting[1][1] += poynting_inst[1][1];
    //    }

    //}//if Ymax

}

void ElectroMagn3DRZ::applyExternalFields(Patch* patch) {
    #ifdef _TODO_RZ            
    #endif
    int imode = 0;

    Field * field;
    for (vector<ExtField>::iterator extfield=extFields.begin(); extfield!=extFields.end(); extfield++ ) {
        string name = LowerCase(extfield->field);
        if      ( El_[imode] && name==LowerCase(El_[imode]->name) ) field = El_[imode];
        else if ( Er_[imode] && name==LowerCase(Er_[imode]->name) ) field = Er_[imode];
        else if ( Et_[imode] && name==LowerCase(Et_[imode]->name) ) field = Et_[imode];
        else if ( Bl_[imode] && name==LowerCase(Bl_[imode]->name) ) field = Bl_[imode];
        else if ( Br_[imode] && name==LowerCase(Br_[imode]->name) ) field = Br_[imode];
        else if ( Bt_[imode] && name==LowerCase(Bt_[imode]->name) ) field = Bt_[imode];
        else field = NULL;

        if( field ) {
            applyExternalField( field, extfield->profile, patch );
        }
    }
    Bl_m[imode]->copyFrom(Bl_[imode]);
    Br_m[imode]->copyFrom(Br_[imode]);
    Bt_m[imode]->copyFrom(Bt_[imode]);
}

void ElectroMagn3DRZ::applyExternalField(Field* my_field,  Profile *profile, Patch* patch) {
    
    cField2D* field2D=static_cast<cField2D*>(my_field);
    
    vector<double> pos(2);
    pos[0]      = dl*((double)(patch->getCellStartingGlobalIndex(0))+(field2D->isDual(0)?-0.5:0.));
    double pos1 = dr*((double)(patch->getCellStartingGlobalIndex(1))+(field2D->isDual(1)?-0.5:0.));
    int N0 = (int)field2D->dims()[0];
    int N1 = (int)field2D->dims()[1];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    for (int i=0 ; i<N0 ; i++) {
        pos[1] = pos1;
        for (int j=0 ; j<N1 ; j++) {
            (*field2D)(i,j) += profile->valueAt(pos);
            pos[1] += dr;
        }
        pos[0] += dl;
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

//! Fold EM fields modes correctly around axis
void ElectroMagn3DRZ::fold_fields(bool diag_flag)
{  

    // Are static casts really necesary here ?

    if (isYmin){

         cField2D* JlRZ ;
         cField2D* JrRZ ;
         cField2D* JtRZ ;

         for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {

             //static cast of the total currents and densities
             JlRZ    = Jl_[imode];
             JrRZ    = Jr_[imode];
             JtRZ    = Jt_[imode];
             for (unsigned int i=0; i<nl_d; i++){
                 for (unsigned int j=0; j<oversize[1]; j++)
                     (*JlRZ)(i,2*oversize[1]-j)+= (*JlRZ)(i,j) ;
             } 
             for (unsigned int i=0; i<nl_p; i++){
                 for (unsigned int j=0; j<oversize[1]; j++)
                     (*JtRZ)(i,2*oversize[1]-j)+= (*JtRZ)(i,j) ;
             }
             if (imode==0){
                 for (unsigned int i=0; i<nl_p; i++){
                      for (unsigned int j=0; j<oversize[1]; j++)
                         (*JrRZ)(i,2*oversize[1]+1-j)+= (*JrRZ)(i,j) ;
                      (*JrRZ)(i,oversize[1]+1)= -(*JrRZ)(i,oversize[1]) ;
                      if (abs((*JrRZ)(i,oversize[1]+1)+ (*JrRZ)(i,oversize[1]))*100!=0.){
                      MESSAGE("careful! "<<abs((*JrRZ)(i,oversize[1]+1)+ (*JrRZ)(i,oversize[1]))*100000 )
                      }
                 }
             }
             else{
                 for (unsigned int i=0; i<nl_p; i++){
                      for (unsigned int j=0; j<oversize[1]+1; j++)
                         (*JrRZ)(i,2*oversize[1]+1-j)+= (*JrRZ)(i,j) ;
                 }
             }
         }
     
         if(diag_flag){
             //Loop on modes for rho
             for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
                 cField2D* rhoRZ   = static_cast<cField2D*>(rho_RZ_[imode]);
                 for (unsigned int i=0; i<nl_p; i++){
                     for (unsigned int j=0; j<oversize[1]; j++)
                         (*rhoRZ)(i,2*oversize[1]-j)+= (*rhoRZ)(i,j) ;
                 } 
             }

             //Loop on all modes and species for J_s
             for (unsigned int ism=0; ism <  n_species*nmodes; ism++){
                 JlRZ    = Jl_s[ism];
                 if ( JlRZ != NULL ) {
                     for (unsigned int i=0; i<nl_d; i++){
                         for (unsigned int j=0; j<oversize[1]; j++){
                             (*JlRZ)(i,2*oversize[1]-j)+= (*JlRZ)(i,j) ;
                         }
                     }
                 }
                 JtRZ    = Jt_s[ism];
                 if ( JtRZ != NULL ) {
                     for (unsigned int i=0; i<nl_p; i++){
                         for (unsigned int j=0; j<oversize[1]; j++)
                             (*JtRZ)(i,2*oversize[1]-j)+= (*JtRZ)(i,j) ;
                     }
                 }
                 JrRZ    = Jr_s[ism];
                 if ( JrRZ != NULL ) {
                     for (unsigned int i=0; i<nl_p; i++){
                         for (unsigned int j=0; j<oversize[1]; j++)
                             (*JrRZ)(i,2*oversize[1]+1-j)+= (*JrRZ)(i,j) ;
                         (*JrRZ)(i,oversize[1]+1)= -(*JrRZ)(i,oversize[1]) ;
                         
                     }
                 }
             }
         }

    } 
    //MESSAGE("folding fields");
    return;
}

//! Evaluating EM fields modes correctly on axis
void ElectroMagn3DRZ::on_axis_fields(bool diag_flag)
{  

    // Are static casts really necesary here ?

    if (isYmin){

         cField2D* JlRZ ;
         cField2D* JrRZ ;
         cField2D* JtRZ ;

         for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {

             //static cast of the total currents and densities
             //JlRZ    = static_cast<cField2D*>(Jl_[imode]);
             //JrRZ    = static_cast<cField2D*>(Jr_[imode]);
             //JtRZ    = static_cast<cField2D*>(Jt_[imode]);
             JlRZ    = Jl_[imode];
             JrRZ    = Jr_[imode];
             JtRZ    = Jt_[imode];
             
             if (imode ==0){
                 for (unsigned int i=0; i<nl_p; i++)
                     (*JtRZ)(i,oversize[1]) = 0.;

                 for (unsigned int i=0; i<nl_p; i++)
                      (*JrRZ)(i,oversize[1])= 0. ;

             }
             else if (imode==1){
                 for (unsigned int i=0; i<nl_p; i++)
                     (*JtRZ)(i,oversize[1]) = - 1./3.* (4.* Icpx * (*JrRZ)(i,oversize[1]+1) + (*JtRZ)(i,oversize[1]));

                 for (unsigned int i=0; i<nl_p; i++)
                      (*JrRZ)(i,oversize[1])= 2.*Icpx* (*JtRZ)(i,oversize[1])-(*JrRZ)(i,oversize[1]+1) ;
             }
         } 
         if(diag_flag){
             for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
                 cField2D* rhoRZ   = static_cast<cField2D*>(rho_RZ_[imode]); 
                 if (imode ==0){
                     for (unsigned int ism=0; ism <  n_species*nmodes; ism++){
                         JtRZ    = Jt_s[ism];
                         if ( JtRZ != NULL ) { 
                             for (unsigned int i=0; i<nl_p; i++){
                                 (*JtRZ)(i,oversize[1]) = 0.;
                             }
                         }
                         JrRZ    = Jr_s[ism];
                         if ( JrRZ != NULL ) { 
                             for (unsigned int i=0; i<nl_p; i++){ 
                                 (*JrRZ)(i,oversize[1])= 0. ;
                             }
                         }
                      }
                 }
                 else if (imode == 1){
                     for (unsigned int i=0; i<nl_p; i++){
                         for (unsigned int j=0; j<oversize[1]; j++)
                             (*rhoRZ)(i,2*oversize[1]-j)+= (*rhoRZ)(i,j) ;
                      }
                     //Loop on all modes and species for J_s
                     for (unsigned int ism=0; ism <  n_species*nmodes; ism++){
                         JtRZ    = Jt_s[ism];
                         if ( JtRZ != NULL ) {
                             for (unsigned int i=0; i<nl_p; i++){
                                 (*JtRZ)(i,oversize[1]) = - 1./3.* (4.* Icpx * (*JrRZ)(i,oversize[1]+1) + (*JtRZ)(i,oversize[1]));
                             }
                         }
                         JrRZ    = Jr_s[ism];
                         if ( JrRZ != NULL ) {
                             for (unsigned int i=0; i<nl_p; i++){
                                 (*JrRZ)(i,oversize[1])= 2.*Icpx* (*JtRZ)(i,oversize[1])-(*JrRZ)(i,oversize[1]+1) ;
                             } 
                         }
                      }
                   }
              }
            }
        } 
    //MESSAGE("folding fields");
    return;
}

