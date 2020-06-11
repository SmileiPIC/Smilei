#include "ElectroMagnAM.h"

#include <cmath>

#include <iostream>
#include <sstream>
#include <complex>
#include "dcomplex.h"
#include "Params.h"
#include "Field2D.h"
#include "cField2D.h"
#include "FieldFactory.h"

#include "Patch.h"
#include "PatchAM.h"
#include <cstring>

#include "Profile.h"

#include "ElectroMagnBC.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for ElectromagnAM
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagnAM::ElectroMagnAM( Params &params, DomainDecomposition *domain_decomposition, vector<Species *> &vecSpecies, Patch *patch ) :
    ElectroMagn( params, domain_decomposition, vecSpecies, patch ),
    isYmin( patch->isYmin() ),
    isYmax( patch->isYmax() )
{

    initElectroMagnAMQuantities( params, patch );
    
    // Charge currents currents and density for each species
    for( unsigned int imode=0; imode<nmodes; imode++ ) {
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
            ostringstream species_mode_name( "" );
            species_mode_name << vecSpecies[ispec]->name_ << "_mode_" << imode;
            Jl_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jl_" , vecSpecies[ispec]->name_, "_mode_",  imode ).c_str(), dimPrim, params );
            Jr_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jr_" , vecSpecies[ispec]->name_, "_mode_",  imode ).c_str(), dimPrim, params );
            Jt_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jt_" , vecSpecies[ispec]->name_, "_mode_",  imode ).c_str(), dimPrim, params );
            rho_AM_s[imode*n_species+ispec] = new cField2D( ( "Rho_"+ species_mode_name.str() ).c_str(), dimPrim );
            if ((imode == 0) && (params.Laser_Envelope_model )){
                Env_Chi_s[ispec] = new Field2D( ( "Env_Chi_"+ species_mode_name.str() ).c_str(), dimPrim );
            }
        }
    }
}//END constructor Electromagn3D


ElectroMagnAM::ElectroMagnAM( ElectroMagnAM *emFields, Params &params, Patch *patch ) :
    ElectroMagn( emFields, params, patch ),
    isYmin( patch->isYmin() ),
    isYmax( patch->isYmax() )
{

    initElectroMagnAMQuantities( params, patch );
    
    // Charge currents currents and density for each species
    for( unsigned int imode=0; imode<nmodes; imode++ ) {
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        
            int ifield = imode*n_species+ispec;

            if( emFields->Jl_s[ifield] != NULL ) {
                if( emFields->Jl_s[ifield]->cdata_ != NULL ) {
                    Jl_s[ifield]  = FieldFactory::createComplex(dimPrim, 0, false, emFields->Jl_s[ifield]->name, params);
                } else {
                    Jl_s[ifield]  = FieldFactory::createComplex(emFields->Jl_s[ifield]->name, dimPrim, params);
                }
            }
            if( emFields->Jr_s[ifield] != NULL ) {
                if( emFields->Jr_s[ifield]->cdata_ != NULL ) {
                    Jr_s[ifield]  = FieldFactory::createComplex(dimPrim, 1, false, emFields->Jr_s[ifield]->name, params);
                } else {
                    Jr_s[ifield]  = FieldFactory::createComplex(emFields->Jr_s[ifield]->name, dimPrim, params);
                }
            }
            if( emFields->Jt_s[ifield] != NULL ) {
                if( emFields->Jt_s[ifield]->cdata_ != NULL ) {
                    Jt_s[ifield]  = FieldFactory::createComplex(dimPrim, 2, false, emFields->Jt_s[ifield]->name, params);
                } else {
                    Jt_s[ifield]  = FieldFactory::createComplex(emFields->Jt_s[ifield]->name, dimPrim, params);
                }
            }
            if( emFields->rho_AM_s[ifield] != NULL ) {
                if( emFields->rho_AM_s[ifield]->cdata_ != NULL ) {
                    rho_AM_s[ifield] = new cField2D( dimPrim, emFields->rho_AM_s[ifield]->name );
                } else {
                    rho_AM_s[ifield]  = new cField2D( emFields->rho_AM_s[ifield]->name, dimPrim );
                }
            }

            if( (imode==0) && (params.Laser_Envelope_model) && (emFields->Env_Chi_s[ifield] != NULL) ) {
                if( emFields->Env_Chi_s[ifield]->data_ != NULL ) {
                    Env_Chi_s[ifield] = new Field2D( dimPrim, emFields->Env_Chi_s[ifield]->name );
                } else {
                    Env_Chi_s[ifield]  = new Field2D( emFields->Env_Chi_s[ifield]->name, dimPrim );
                }    
            }

        }
        
    }
    
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn3D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::initElectroMagnAMQuantities( Params &params, Patch *patch )
{

    nmodes = params.nmodes;

    PatchAM *patchAM = static_cast<PatchAM *>( patch );
    invR = &(patchAM->invR[0]);
    invRd = &(patchAM->invRd[0]);
    
    // Species charge currents and density
    Jl_s.resize( n_species*nmodes );
    Jr_s.resize( n_species*nmodes );
    Jt_s.resize( n_species*nmodes );
    rho_AM_s.resize( n_species*nmodes );
    for( unsigned int ispec=0; ispec<n_species*nmodes; ispec++ ) {
        Jl_s[ispec]  = NULL;
        Jr_s[ispec]  = NULL;
        Jt_s[ispec]  = NULL;
        rho_AM_s[ispec] = NULL;
    }

    if (params.Laser_Envelope_model){
        Env_Chi_s.resize( n_species );
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        Env_Chi_s[ispec]  = NULL;
        }
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
    j_glob_ = patch->getCellStartingGlobalIndex( 1 );
    
    
    // ----------------------
    // Electromagnetic fields
    // ----------------------
    
    dimPrim.resize( nDim_field );
    dimDual.resize( nDim_field );
    
    // Dimension of the primal and dual grids
    for( size_t i=0 ; i<nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2-(params.is_pxr);
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nl_p = n_space[0]+1+2*oversize[0];
    nl_d = n_space[0]+2+2*oversize[0]-(params.is_pxr);
    // number of nodes of the primal and dual grid in the y-direction
    nr_p = n_space[1]+1+2*oversize[1];
    nr_d = n_space[1]+2+2*oversize[1]-(params.is_pxr);
    
    // Allocation of the EM fields
    
    El_.resize( nmodes );
    Er_.resize( nmodes );
    Et_.resize( nmodes );
    Bl_.resize( nmodes );
    Br_.resize( nmodes );
    Bt_.resize( nmodes );
    Bl_m.resize( nmodes );
    Br_m.resize( nmodes );
    Bt_m.resize( nmodes );
    
    // Total charge currents and densities
    Jl_.resize( nmodes );
    Jr_.resize( nmodes );
    Jt_.resize( nmodes );
    rho_AM_.resize( nmodes );

    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        ostringstream mode_id( "" );
        mode_id << "_mode_" << imode;
        
        El_[imode]  = FieldFactory::createComplex( dimPrim, 0, false, ( "El"+mode_id.str() ).c_str(), params );
        Er_[imode]  = FieldFactory::createComplex( dimPrim, 1, false, ( "Er"+mode_id.str() ).c_str(), params );
        Et_[imode]  = FieldFactory::createComplex( dimPrim, 2, false, ( "Et"+mode_id.str() ).c_str(), params );
        Bl_[imode]  = FieldFactory::createComplex( dimPrim, 0, true, ( "Bl"+mode_id.str() ).c_str(), params );
        Br_[imode]  = FieldFactory::createComplex( dimPrim, 1, true, ( "Br"+mode_id.str() ).c_str(), params );
        Bt_[imode]  = FieldFactory::createComplex( dimPrim, 2, true, ( "Bt"+mode_id.str() ).c_str(), params );
        Bl_m[imode] = FieldFactory::createComplex( dimPrim, 0, true, ( "Bl_m"+mode_id.str() ).c_str(), params );
        Br_m[imode] = FieldFactory::createComplex( dimPrim, 1, true, ( "Br_m"+mode_id.str() ).c_str(), params );
        Bt_m[imode] = FieldFactory::createComplex( dimPrim, 2, true, ( "Bt_m"+mode_id.str() ).c_str(), params );
        
        // Total charge currents and densities
        Jl_[imode]   = FieldFactory::createComplex( dimPrim, 0, false, ( "Jl"+mode_id.str() ).c_str(), params );
        Jr_[imode]   = FieldFactory::createComplex( dimPrim, 1, false, ( "Jr"+mode_id.str() ).c_str(), params );
        Jt_[imode]   = FieldFactory::createComplex( dimPrim, 2, false, ( "Jt"+mode_id.str() ).c_str(), params );
        rho_AM_[imode]  = new cField2D( dimPrim, ( "Rho"+mode_id.str() ).c_str() );
    }

    if(params.is_pxr == true) {
        rho_old_AM_.resize( nmodes );
        for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
            ostringstream mode_id( "" );
            mode_id << "_mode_" << imode;
            rho_old_AM_[imode]  = new cField2D( dimPrim, ( "RhoOld"+mode_id.str() ).c_str() );
        }
    } 

    

    if( params.Laser_Envelope_model ) {
        Env_A_abs_  = new Field2D( dimPrim, "Env_A_abs_mode_0" );
        Env_Chi_    = new Field2D( dimPrim, "Env_Chi_mode_0" );
        Env_E_abs_  = new Field2D( dimPrim, "Env_E_abs_mode_0" );
        Env_Ex_abs_ = new Field2D( dimPrim, "Env_Ex_abs_mode_0" );
    }
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for( unsigned int i=0 ; i<nDim_field ; i++ ) {
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
    for( unsigned int i=0 ; i<3 ; i++ )
        for( unsigned int isDual=0 ; isDual<2 ; isDual++ ) {
            istart[i][isDual] = 0;
        }
    for( unsigned int i=0 ; i<nDim_field ; i++ ) {
        for( unsigned int isDual=0 ; isDual<2 ; isDual++ ) {
            istart[i][isDual] = oversize[i];
            if( patch->Pcoordinates[i]!=0 ) {
                istart[i][isDual]+=1;
            }
        }
    }
    
    // bufsize = nelements
    for( unsigned int i=0 ; i<3 ; i++ )
        for( unsigned int isDual=0 ; isDual<2 ; isDual++ ) {
            bufsize[i][isDual] = 1;
        }
        
    for( unsigned int i=0 ; i<nDim_field ; i++ ) {
        for( int isDual=0 ; isDual<2 ; isDual++ ) {
            bufsize[i][isDual] = n_space[i] + 1;
        }
        
        for( int isDual=0 ; isDual<2 ; isDual++ ) {
            bufsize[i][isDual] += isDual;
            if( params.number_of_patches[i]!=1 ) {
            
                if( ( !isDual ) && ( patch->Pcoordinates[i]!=0 ) ) {
                    bufsize[i][isDual]--;
                } else if( isDual ) {
                    bufsize[i][isDual]--;
                    if( ( patch->Pcoordinates[i]!=0 ) && ( patch->Pcoordinates[i]!=params.number_of_patches[i]-1 ) ) {
                        bufsize[i][isDual]--;
                    }
                }
                
            } // if ( params.number_of_patches[i]!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field
}


void ElectroMagnAM::finishInitialization( int nspecies, Patch *patch )
{
    // Fill allfields
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
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
        allFields.push_back( rho_AM_[imode] );
        if( (imode ==0) && (Env_A_abs_ != NULL) ) {
            allFields.push_back( Env_A_abs_ );
            allFields.push_back( Env_Chi_ );
            allFields.push_back( Env_E_abs_ );
            allFields.push_back( Env_Ex_abs_ );
        }
    }
    
    // For species-related fields
    // The order is necessary in DiagnosticProbes - DO NOT CHANGE -
    species_starts.resize( 0 );
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        species_starts.push_back( allFields.size() );
        for( unsigned int imode=0; imode<nmodes; imode++ ) {
            int ifield = imode*n_species+ispec;
            allFields.push_back( Jl_s[ifield] );
            allFields.push_back( Jr_s[ifield] );
            allFields.push_back( Jt_s[ifield] );
            allFields.push_back( rho_AM_s[ifield] );
        }
        if( Env_A_abs_ != NULL ){
            int ifield = 0*n_species+ispec; // only mode 0
            allFields.push_back( Env_Chi_s[ifield] );
        }
    }
    species_starts.push_back( allFields.size() );
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for ElectromagnAM
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagnAM::~ElectroMagnAM()
{
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        delete El_[imode];
        delete Er_[imode];
        delete Et_[imode];
        delete Bl_[imode];
        delete Br_[imode];
        delete Bt_[imode];
        if (!is_pxr) {
            delete Bl_m[imode];
            delete Br_m[imode];
            delete Bt_m[imode];
        }
        
        delete Jl_[imode];
        delete Jr_[imode];
        delete Jt_[imode];
        delete rho_AM_[imode];
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
            int ifield = imode*n_species+ispec;
            delete Jl_s[ifield];
            delete Jr_s[ifield];
            delete Jt_s[ifield];
            delete rho_AM_s[ifield];            
        }
     
    }
    
}//END ElectroMagnAM


void ElectroMagnAM::restartRhoJ()
{
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        Jl_[imode] ->put_to( 0. );
        Jr_[imode] ->put_to( 0. );
        Jt_[imode] ->put_to( 0. );
        rho_AM_[imode]->put_to( 0. );
    }
}

void ElectroMagnAM::restartRhoJs()
{
    for( unsigned int ispec=0 ; ispec < n_species*nmodes ; ispec++ ) {
        if( Jl_s [ispec] ) {
            Jl_s [ispec]->put_to( 0. );
        }
        if( Jr_s [ispec] ) {
            Jr_s [ispec]->put_to( 0. );
        }
        if( Jt_s [ispec] ) {
            Jt_s [ispec]->put_to( 0. );
        }
        if( rho_AM_s[ispec] ) {
            rho_AM_s[ispec]->put_to( 0. );
        }
    }
    
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        Jl_[imode] ->put_to( 0. );
        Jr_[imode]->put_to( 0. );
        Jt_[imode]->put_to( 0. );
        rho_AM_[imode]->put_to( 0. );
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


void ElectroMagnAM::initPoisson( Patch *patch )
{

    
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------
    
    index_min_p_.resize( 2, 0 );
    index_max_p_.resize( 2, 0 );
    
    index_min_p_[0] = oversize[0];
    index_min_p_[1] = oversize[1];
    index_max_p_[0] = nl_p - 2 - oversize[0];
    index_max_p_[1] = nr_p - 2 - oversize[1];
    if( patch->isXmin() ) {
        index_min_p_[0] = 0;
    }
    if( patch->isXmax() ) {
        index_max_p_[0] = nl_p-1;
    }
    
    phi_AM_ = new cField2D( dimPrim );  // scalar potential
    r_AM_   = new cField2D( dimPrim );  // residual vector
    p_AM_   = new cField2D( dimPrim );  // direction vector
    Ap_AM_  = new cField2D( dimPrim );  // A*p vector
    
} // initPoisson

void ElectroMagnAM::initPoisson_init_phi_r_p_Ap( Patch *patch, unsigned int imode ){

    cField2D *rho   = rho_AM_[imode];
    double j_;
    double dr_sq_dl = dr*dr*dl; // Poisson's equation including the source term is multiplied by r_j*dr*dl
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            j_ = (double)( j_glob_+j);
            ( *phi_AM_ )( i, j )   = 0.; 
            ( *r_AM_ )( i, j )     = -(( *rho )( i, j ))*j_*dr_sq_dl; 
            ( *p_AM_ )( i, j )     = ( *r_AM_ )( i, j );
        }//j
    }//i

}


double ElectroMagnAM::compute_r()
{
    double rnew_dot_rnew_localAM_( 0. );
    for( unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++ ) {
        for( unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++ ) {
            rnew_dot_rnew_localAM_ += std::abs(( *r_AM_ )( i, j ))*std::abs(( *r_AM_ )( i, j ));
        }
    }
    return rnew_dot_rnew_localAM_;
} // compute_r

void ElectroMagnAM::compute_Ap( Patch *patch )
{
#ifdef _TODO_AM
#endif
} // compute_pAp

void ElectroMagnAM::compute_Ap_relativistic_Poisson_AM( Patch *patch, double gamma_mean, unsigned int imode )
{
    
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
      
    // relativistic Poisson's equation in finite differences is multiplied by r_j*dr*dl to condition it before conjugate gradient

    double dr_sq_ov_dl_ov_gamma_sq   = ( dr*dr )/dl/( gamma_mean*gamma_mean );
//    double dl_ov_2                   = dl/2.;
    double m_sq_dl                   = (double)(imode*imode)*dl;
    double j_;
    unsigned int j_min =max(1,isYmin*3); // prevent a segmentation fault
    unsigned int i_min =1; 
    unsigned int i_max = nl_p-1; 
    unsigned int j_max = nr_p-1; 
   
    // vector product Ap = A*p
    for( unsigned int i=i_min; i<i_max; i++ ) {
        for( unsigned int j=j_min; j<j_max; j++ ) {
            j_ = (double)( j_glob_+j);
            ( *Ap_AM_ )( i, j )= j_ * dr_sq_ov_dl_ov_gamma_sq * (          ( *p_AM_ )( i-1, j   )-2.*   ( *p_AM_ )( i, j   )+         ( *p_AM_ )( i+1, j ) )
                               + dl                           * ( (j_-0.5)*( *p_AM_ )( i  , j-1 )-2.*j_*( *p_AM_ )( i, j   )+(j_+0.5)*( *p_AM_ )( i, j+1 ) )
                               - m_sq_dl/j_                  *                                          ( *p_AM_ )( i, j   );                     
        }//j
    }//i
    
    
    // Axis BC
    if( patch->isYmin() ) {
        unsigned int j=2;
        j_ = (double)( j_glob_+j+0.5);
        for( unsigned int i=i_min; i<i_max; i++ ) { // radial and azimuthal derivative are zero on axis r=0 (p = phi is all on primal grid)
            ( *Ap_AM_ )( i, j )= j_ * dr_sq_ov_dl_ov_gamma_sq * (          ( *p_AM_ )( i-1, j   )-2.*   ( *p_AM_ )( i, j   )+         ( *p_AM_ )( i+1, j ) )
                               + j_ * dl * 2.                 * (                                       ( *p_AM_ )( i, j+1 )-         ( *p_AM_ )( i  , j)  );                           
        }
    }

    // Xmin BC
    if( patch->isXmin() ) { // p = phi = 0 on the left border
        for( unsigned int j=1; j<j_max; j++ ) {
            ( *Ap_AM_ )( 0, j )     = 0.;
        }
        // at corners
        ( *Ap_AM_ )( 0, 0 )          = 0.;
        ( *Ap_AM_ )( 0, nr_p-1 )     = 0.;
    }
    
    // Xmax BC
    if( patch->isXmax() ) { // p = phi = 0 on the right border 
    
        for( unsigned int j=1; j<j_max; j++ ) {
            ( *Ap_AM_ )( nl_p-1, j )= 0.;
        }
        // at corners
        ( *Ap_AM_ )( nl_p-1, 0 )     = 0.;
        ( *Ap_AM_ )( nl_p-1, nr_p-1 )= 0.;
    }
    

} // compute_pAp

void ElectroMagnAM::compute_Ap_Poisson_AM( Patch *patch, unsigned int imode )
{
      
    // Poisson's equation in finite differences is multiplied by r_j*dr*dl to condition it before conjugate gradient

    double dr_sq_ov_dl   = ( dr*dr )/dl;
//    double dl_ov_2                   = dl/2.;
    double m_sq_dl       = (double)(imode*imode)*dl;
    double j_;
    unsigned int j_min =max(1,isYmin*3); // prevent a segmentation fault
    unsigned int i_min =1; 
    unsigned int i_max = nl_p-1; 
    unsigned int j_max = nr_p-1; 
   
    // vector product Ap = A*p
    for( unsigned int i=i_min; i<i_max; i++ ) {
        for( unsigned int j=j_min; j<j_max; j++ ) {
            j_ = (double)( j_glob_+j);
            ( *Ap_AM_ )( i, j )= j_ * dr_sq_ov_dl             * (          ( *p_AM_ )( i-1, j   )-2.*   ( *p_AM_ )( i, j   )+         ( *p_AM_ )( i+1, j ) )
                               + dl                           * ( (j_-0.5)*( *p_AM_ )( i  , j-1 )-2.*j_*( *p_AM_ )( i, j   )+(j_+0.5)*( *p_AM_ )( i, j+1 ) )
                               - m_sq_dl/j_                  *                                          ( *p_AM_ )( i, j   );                     
        }//j
    }//i
    
    
    // Axis BC
    if( patch->isYmin() ) {
        unsigned int j=2;
        j_ = (double)( j_glob_+j+0.5);
        for( unsigned int i=i_min; i<i_max; i++ ) { // radial and azimuthal derivative are zero on axis r=0 (p = phi is all on primal grid)
            ( *Ap_AM_ )( i, j )= j_ * dr_sq_ov_dl             * (          ( *p_AM_ )( i-1, j   )-2.*   ( *p_AM_ )( i, j   )+         ( *p_AM_ )( i+1, j ) )
                               + j_ * dl * 2.                 * (                                       ( *p_AM_ )( i, j+1 )-         ( *p_AM_ )( i  , j)  );                           
        }
    }

    // Xmin BC
    if( patch->isXmin() ) { // p = phi = 0 on the left border
        for( unsigned int j=1; j<j_max; j++ ) {
            ( *Ap_AM_ )( 0, j )     = 0.;
        }
        // at corners
        ( *Ap_AM_ )( 0, 0 )          = 0.;
        ( *Ap_AM_ )( 0, nr_p-1 )     = 0.;
    }
    
    // Xmax BC
    if( patch->isXmax() ) { // p = phi = 0 on the right border 
    
        for( unsigned int j=1; j<j_max; j++ ) {
            ( *Ap_AM_ )( nl_p-1, j )= 0.;
        }
        // at corners
        ( *Ap_AM_ )( nl_p-1, 0 )     = 0.;
        ( *Ap_AM_ )( nl_p-1, nr_p-1 )= 0.;
    }
    

} // compute_pAp

std::complex<double> ElectroMagnAM::compute_pAp_AM()
{
    std::complex<double> p_dot_Ap_local = 0.0;
    for( unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++ ) {
        for( unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++ ) {
            p_dot_Ap_local += ( *p_AM_ )( i, j )*std::conj(( *Ap_AM_ )( i, j ));
        }
    }
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagnAM::update_pand_r_AM( double r_dot_r, std::complex<double> p_dot_Ap )
{
    std::complex<double> alpha_k = r_dot_r/p_dot_Ap;
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *phi_AM_ )( i, j ) += alpha_k * ( *p_AM_ )( i, j );
            ( *r_AM_ )( i, j )   -= alpha_k * ( *Ap_AM_ )( i, j );
        }
    }
    
} // update_pand_r

void ElectroMagnAM::update_p( double rnew_dot_rnew, double r_dot_r )
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *p_AM_ )( i, j ) = ( *r_AM_ )( i, j ) + beta_k * ( *p_AM_ )( i, j );
        }
    }
} // update_p



void ElectroMagnAM::initRelativisticPoissonFields( Patch *patch ){
    // ------ Init temporary fields for relativistic field initialization
    
    // E fields centered as in FDTD, to be added to the already present electric fields
    El_rel_  = new cField2D( dimPrim, 0, false, "El_rel" );
    Er_rel_  = new cField2D( dimPrim, 1, false, "Er_rel" );
    Et_rel_  = new cField2D( dimPrim, 2, false, "Et_rel" );
    
    
    // B fields centered as the E fields in FDTD (Bx null)
    Bl_rel_  = new cField2D( dimPrim, 0, true,  "Bl_rel" );  // null
    Br_rel_  = new cField2D( dimPrim, 2, false,  "Br_rel" ); // centered as Et initially
    Bt_rel_  = new cField2D( dimPrim, 1, false,  "Bt_rel" ); // centered as Er initially
    
    
    // ----- B fields centered as in FDTD, to be added to the already present magnetic fields
    
    // B field advanced by dt/2
    Bl_rel_t_plus_halfdt_  = new cField2D( dimPrim, 0, true,  "Bl_rel_t_plus_halfdt" );
    Br_rel_t_plus_halfdt_  = new cField2D( dimPrim, 1, true,  "Br_rel_t_plus_halfdt" );
    Bt_rel_t_plus_halfdt_  = new cField2D( dimPrim, 2, true,  "Bt_rel_t_plus_halfdt" );
    // B field "advanced" by -dt/2
    Bl_rel_t_minus_halfdt_  = new cField2D( dimPrim, 0, true,  "Bl_rel_t_plus_halfdt" );
    Br_rel_t_minus_halfdt_  = new cField2D( dimPrim, 1, true,  "Br_rel_t_plus_halfdt" );
    Bt_rel_t_minus_halfdt_  = new cField2D( dimPrim, 2, true,  "Bt_rel_t_plus_halfdt" );


}

void ElectroMagnAM::initPoissonFields( Patch *patch ){
    // ------ Init temporary fields for field initialization
    
    // E fields centered as in FDTD, to be added to the already present electric fields
    El_Poisson_  = new cField2D( dimPrim, 0, false, "El_Poisson" );
    Er_Poisson_  = new cField2D( dimPrim, 1, false, "Er_Poisson" );
    Et_Poisson_  = new cField2D( dimPrim, 2, false, "Et_Poisson" );

}

void ElectroMagnAM::initE_relativistic_Poisson_AM( Patch *patch, double gamma_mean, unsigned int imode )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    cField2D *ElAM  = static_cast<cField2D *>( El_rel_ );
    cField2D *ErAM  = static_cast<cField2D *>( Er_rel_ );
    cField2D *EtAM  = static_cast<cField2D *>( Et_rel_ );
    //cField2D *rhoAM = static_cast<cField2D *>( rho_AM_[imode] );

    complex<double>     i1 = std::complex<double>( 0., 1 );

    // ------------------------------------------
    // Compute the fields El, Er and Et
    // ------------------------------------------
    
    // El
    MESSAGE( 1, "Computing El from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_d-1; i++ ) {
        for( unsigned int j=1; j<nr_p; j++ ) {
            ( *ElAM )( i, j ) = ( ( *phi_AM_ )( i-1, j )-( *phi_AM_ )( i, j ) )/dl/gamma_mean/gamma_mean;
        }
    }
    MESSAGE( 1, "El: done" );
    // Er
    MESSAGE( 1, "Computing Er from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_d; j++ ) {
            ( *ErAM )( i, j ) = ( ( *phi_AM_ )( i, j-1 )-( *phi_AM_ )( i, j ) )/dr;
        }
    }
    MESSAGE( 1, "Er: done" );
    // Et
    MESSAGE( 1, "Computing Er from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_p; j++ ) {
            ( *EtAM )( i, j ) = i1*((double )imode)/(((double)( j_glob_+j ))*dr)* ( *phi_AM_ )( i, j );
        }
    }
    MESSAGE( 1, "Et: done" );

    
    if( isYmin ) { // Conditions on axis
        unsigned int j=2;
        if( imode==0 ) {
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *EtAM )( i, j )=0;
            }
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *ErAM )( i, j )= -( *ErAM )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_d ; i++ ) {
                ( *ElAM )( i, j ) = ( *ElAM )( i, j+1 ) ; //( *ElAM )( i, j ) = ( *ElAM )( i, j+1 ) ;  // not sure about this one
            }
        } else if( imode==1 ) {
            for( unsigned int i=0 ; i<nl_d  ; i++ ) {
                ( *ElAM )( i, j )= 0;
            }
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *EtAM )( i, j )= -1./3.*( 4.*i1*( *ErAM )( i, j+1 )+( *EtAM )( i, j+1 ) );
            }
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                ( *ErAM )( i, j )=2.*i1*( *EtAM )( i, j )-( *ErAM )( i, j+1 );
            }
        } else {
            for( unsigned int  i=0 ; i<nl_d; i++ ) {
                ( *ElAM )( i, j )= 0;
            }
            for( unsigned int  i=0 ; i<nl_p; i++ ) {
                ( *ErAM )( i, j )= -( *ErAM )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_p; i++ ) {
                ( *EtAM )( i, j )= 0;
            }
        }
    }
  
    
} // initE_relativistic_Poisson_AM

void ElectroMagnAM::initE_Poisson_AM( Patch *patch, unsigned int imode )
{
    
    cField2D *ElAM  = static_cast<cField2D *>( El_Poisson_ );
    cField2D *ErAM  = static_cast<cField2D *>( Er_Poisson_ );
    cField2D *EtAM  = static_cast<cField2D *>( Et_Poisson_ );

    complex<double>     i1 = std::complex<double>( 0., 1 );

    // ------------------------------------------
    // Compute the fields El, Er and Et
    // ------------------------------------------
    
    // El
    MESSAGE( 1, "Computing El from scalar potential, Poisson problem" );
    for( unsigned int i=1; i<nl_d-1; i++ ) {
        for( unsigned int j=1; j<nr_p; j++ ) {
            ( *ElAM )( i, j ) = ( ( *phi_AM_ )( i-1, j )-( *phi_AM_ )( i, j ) )/dl;
        }
    }
    MESSAGE( 1, "El: done" );
    // Er
    MESSAGE( 1, "Computing Er from scalar potential, Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_d; j++ ) {
            ( *ErAM )( i, j ) = ( ( *phi_AM_ )( i, j-1 )-( *phi_AM_ )( i, j ) )/dr;
        }
    }
    MESSAGE( 1, "Er: done" );
    // Et
    MESSAGE( 1, "Computing Er from scalar potential, Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_p; j++ ) {
            ( *EtAM )( i, j ) = i1*((double )imode)/(((double)( j_glob_+j ))*dr)* ( *phi_AM_ )( i, j );
        }
    }
    MESSAGE( 1, "Et: done" );

    
    if( isYmin ) { // Conditions on axis
        unsigned int j=2;
        if( imode==0 ) {
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *EtAM )( i, j )=0;
            }
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *ErAM )( i, j )= -( *ErAM )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_d ; i++ ) {
                ( *ElAM )( i, j ) = ( *ElAM )( i, j+1 ) ; //( *ElAM )( i, j ) = ( *ElAM )( i, j+1 ) ;  // not sure about this one
            }
        } else if( imode==1 ) {
            for( unsigned int i=0 ; i<nl_d  ; i++ ) {
                ( *ElAM )( i, j )= 0;
            }
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *EtAM )( i, j )= -1./3.*( 4.*i1*( *ErAM )( i, j+1 )+( *EtAM )( i, j+1 ) );
            }
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                ( *ErAM )( i, j )=2.*i1*( *EtAM )( i, j )-( *ErAM )( i, j+1 );
            }
        } else {
            for( unsigned int  i=0 ; i<nl_d; i++ ) {
                ( *ElAM )( i, j )= 0;
            }
            for( unsigned int  i=0 ; i<nl_p; i++ ) {
                ( *ErAM )( i, j )= -( *ErAM )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_p; i++ ) {
                ( *EtAM )( i, j )= 0;
            }
        }
    }
  
    
} // initE_Poisson_AM


void ElectroMagnAM::initB_relativistic_Poisson_AM( Patch *patch, double gamma_mean )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    
    cField2D *ErAM  = static_cast<cField2D *>( Er_rel_ );
    cField2D *EtAM  = static_cast<cField2D *>( Et_rel_ );
    
    cField2D *BlAM  = static_cast<cField2D *>( Bl_rel_ ); // Bl is zero everywhere
    cField2D *BrAM  = static_cast<cField2D *>( Br_rel_ );
    cField2D *BtAM  = static_cast<cField2D *>( Bt_rel_ );

    // ------------------------------------------
    // Compute the field Bl, Br, Bt
    // ------------------------------------------
    
    double beta_mean = sqrt( 1.-1./gamma_mean/gamma_mean );
    MESSAGE( 0, "In relativistic Poisson solver, gamma_mean = " << gamma_mean );
    
    // Bl^(p,d) is identically zero
    MESSAGE( 1, "Computing Bl, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_d-1; j++ ) {
            ( *BlAM )( i, j ) = 0.;
        }
    }
    MESSAGE( 1, "Bl: done" );
    
    // Br^(d,d) from Et^(p,p)
    MESSAGE( 1, "Computing Br from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_p-1; j++ ) {
            ( *BrAM )( i, j ) = -beta_mean*( *EtAM )( i, j );
        }
    }
    MESSAGE( 1, "Br: done" );

    // Bt^(d,p) from Er^(p,d)
    MESSAGE( 1, "Computing Bt from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nl_p-1; i++ ) {
        for( unsigned int j=1; j<nr_d-1; j++ ) {
            ( *BtAM )( i, j ) = beta_mean*( *ErAM )( i, j );
        }
    }
    MESSAGE( 1, "Bt: done" );


    // Should we add BCs here?
    
    
} // initB_relativistic_Poisson_AM

void ElectroMagnAM::center_fields_from_relativistic_Poisson_AM( Patch *patch )
{

    // B field centered in time as E field, at time t
    cField2D *BlAMrel  = static_cast<cField2D *>( Bl_rel_ );
    cField2D *BrAMrel  = static_cast<cField2D *>( Br_rel_ );
    cField2D *BtAMrel  = static_cast<cField2D *>( Bt_rel_ );
    
    // B field centered in time at time t+dt/2
    cField2D *BlAM  = static_cast<cField2D *>( Bl_rel_t_plus_halfdt_ );
    cField2D *BrAM  = static_cast<cField2D *>( Br_rel_t_plus_halfdt_ );
    cField2D *BtAM  = static_cast<cField2D *>( Bt_rel_t_plus_halfdt_ );
    // B field centered in time at time t-dt/2
    cField2D *BlAM0  = static_cast<cField2D *>( Bl_rel_t_minus_halfdt_ );
    cField2D *BrAM0  = static_cast<cField2D *>( Br_rel_t_minus_halfdt_ );
    cField2D *BtAM0  = static_cast<cField2D *>( Bt_rel_t_minus_halfdt_ );
    
    
    // The B_rel fields, centered as B, will be advanced by dt/2 and -dt/2
    // for proper centering in FDTD, but first they have to be centered in space
    // The advance by dt and -dt and the sum to the existing grid fields is performed in
    // ElectroMagn2D::sum_rel_fields_to_em_fields
    
    // Bl (p,d)   Bl_rel is identically zero and centered as Bl, no special interpolation of indices
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_d; j++ ) {
            ( *BlAM )( i, j ) = ( *BlAMrel )( i, j );
            ( *BlAM0 )( i, j )= ( *BlAMrel )( i, j );
        }
    }
    
    // ---------- center the B fields
    // Br (d,p) - remember that Byrel is centered as Etrel (p,p)
    for( unsigned int i=1; i<nl_d-1; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *BrAM )( i, j ) = 0.5 * ( ( *BrAMrel )( i, j ) + ( *BrAMrel )( i-1, j ) );
            ( *BrAM0 )( i, j )= 0.5 * ( ( *BrAMrel )( i, j ) + ( *BrAMrel )( i-1, j ) );
        }
    }
    
    // Bt (d,d) - remember that Btrel is centered as Errel (p,d)
    for( unsigned int i=1; i<nl_d-1; i++ ) {
        for( unsigned int j=0; j<nr_d; j++ ) {
            ( *BtAM )( i, j ) = 0.5 * ( ( *BtAMrel )( i, j ) + ( *BtAMrel )( i-1, j ) );
            ( *BtAM0 )( i, j )= 0.5 * ( ( *BtAMrel )( i, j ) + ( *BtAMrel )( i-1, j ) );
        }
    }
    
}

void ElectroMagnAM::sum_rel_fields_to_em_fields_AM( Patch *patch, Params &params, unsigned int imode )
{
    cField2D *ElAMrel  = static_cast<cField2D *>( El_rel_ );
    cField2D *ErAMrel  = static_cast<cField2D *>( Er_rel_ );
    cField2D *EtAMrel  = static_cast<cField2D *>( Et_rel_ );
    
    // B_t_plus_halfdt
    cField2D *Bl_rel_t_plus_halfdt = static_cast<cField2D *>( Bl_rel_t_plus_halfdt_ );
    cField2D *Br_rel_t_plus_halfdt = static_cast<cField2D *>( Br_rel_t_plus_halfdt_ );
    cField2D *Bt_rel_t_plus_halfdt = static_cast<cField2D *>( Bt_rel_t_plus_halfdt_ );
    
    // B_t_minus_halfdt
    cField2D *Bl_rel_t_minus_halfdt = static_cast<cField2D *>( Bl_rel_t_minus_halfdt_ );
    cField2D *Br_rel_t_minus_halfdt = static_cast<cField2D *>( Br_rel_t_minus_halfdt_ );
    cField2D *Bt_rel_t_minus_halfdt = static_cast<cField2D *>( Bt_rel_t_minus_halfdt_ );
    
    // E and B fields already existing on the grid
    cField2D *ElAM  = static_cast<cField2D *>( El_[imode] );
    cField2D *ErAM  = static_cast<cField2D *>( Er_[imode] );
    cField2D *EtAM  = static_cast<cField2D *>( Et_[imode] );
    cField2D *BlAM  = static_cast<cField2D *>( Bl_[imode] );
    cField2D *BrAM  = static_cast<cField2D *>( Br_[imode]);
    cField2D *BtAM  = static_cast<cField2D *>( Bt_[imode] );
    cField2D *BlAM0  = static_cast<cField2D *>( Bl_m[imode] );
    cField2D *BrAM0  = static_cast<cField2D *>( Br_m[imode] );
    cField2D *BtAM0  = static_cast<cField2D *>( Bt_m[imode] );

    complex<double>     i1 = std::complex<double>( 0., 1 );
    double dt = params.timestep;
    // El (d,p)
    for( unsigned int i=0; i<nl_d; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *ElAM )( i, j ) = ( *ElAM )( i, j ) + ( *ElAMrel )( i, j );
        }
    }
    
    // Er (p,d)
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_d; j++ ) {
            ( *ErAM )( i, j ) = ( *ErAM )( i, j ) + ( *ErAMrel )( i, j );
        }
    }
    
    // Et (p,p)
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *EtAM )( i, j ) = ( *EtAM )( i, j ) + ( *EtAMrel )( i, j );
        }
    }
    
    
    
    // Since Brel is centered in time as E, it is inconsistent with FDTD,
    // where E and B are staggered in time.
    // Possible solution:
    // Use FDTD scheme to integrate Maxwell-Faraday equation forward in time by dt/2 to obtain B
    // Use FDTD scheme to integrate Maxwell-Faraday equation backwards in time by dt/2 to obtain Bm
    // Add the forward-evolved and backward-evolved fields to the grid fields
    
    
    
    // Magnetic field Bl^(p,d)
    for( unsigned int i=0 ; i<nl_p;  i++ ) {
        for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
            // forward advance by dt/2
            ( *Bl_rel_t_plus_halfdt )( i, j ) += - dt/2./( ( j_glob_+j-0.5 )*dr ) * ( ( double )( j+j_glob_ )*( *EtAMrel )( i, j ) 
                                                 - ( double )( j+j_glob_-1. )*( *EtAMrel )( i, j-1 ) + i1*( double )imode*( *ErAMrel )( i, j ) );
            
            // backward advance by dt/2
            ( *Bl_rel_t_minus_halfdt )( i, j )-= - dt/2./( ( j_glob_+j-0.5 )*dr ) * ( ( double )( j+j_glob_ )*( *EtAMrel )( i, j ) 
                                                 - ( double )( j+j_glob_-1. )*( *EtAMrel )( i, j-1 ) + i1*( double )imode*( *ErAMrel )( i, j ) );
        }
    }
    
    // Magnetic field Br^(d,p)
    for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
        for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
            // forward advance by dt/2
            ( *Br_rel_t_plus_halfdt )( i, j ) += dt_ov_dl/2. * ( ( *EtAMrel )( i, j ) - ( *EtAMrel )( i-1, j ) )
                                                 +i1*dt/2.*( double )imode/( ( double )( j_glob_+j )*dr )*( *ElAMrel )( i, j ) ;
            // backward advance by dt/2
            ( *Br_rel_t_minus_halfdt )( i, j )-= dt_ov_dl/2. * ( ( *EtAMrel )( i, j ) - ( *EtAMrel )( i-1, j ) )
                                                 +i1*dt/2.*( double )imode/( ( double )( j_glob_+j )*dr )*( *ElAMrel )( i, j ) ;
        }
    }
    
    // Magnetic field Bt^(d,d)
    for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
        for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
            // forward advance by dt/2
            ( *Bt_rel_t_plus_halfdt )( i, j ) += dt_ov_dr/2. * ( ( *ElAMrel )( i, j ) - ( *ElAMrel )( i, j-1 ) )
                                                -dt_ov_dl/2. * ( ( *ErAMrel )( i, j ) - ( *ErAMrel )( i-1, j ) );
            
            // backward advance by dt/2
            ( *Bt_rel_t_minus_halfdt )( i, j )-= dt_ov_dr/2. * ( ( *ElAMrel )( i, j ) - ( *ElAMrel )( i, j-1 ) )
                                                -dt_ov_dl/2. * ( ( *ErAMrel )( i, j ) - ( *ErAMrel )( i-1, j ) );
        }
    }
    

    // Boundary conditions on Axis
    if( isYmin ) {
        unsigned int j=2;
        if( imode==0 ) {
            for( unsigned int i=0 ; i<nl_d ; i++ ) {
                ( *Br_rel_t_plus_halfdt )( i, j ) =0;
                ( *Br_rel_t_minus_halfdt )( i, j )=0;
            }
            for( unsigned int i=0 ; i<nl_d ; i++ ) {
                ( *Bt_rel_t_plus_halfdt )( i, j ) = -( *Bt_rel_t_plus_halfdt )( i, j+1 );
                ( *Bt_rel_t_minus_halfdt )( i, j )= -( *Bt_rel_t_minus_halfdt )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                ( *Bl_rel_t_plus_halfdt )( i, j ) = ( *Bl_rel_t_plus_halfdt )( i, j+1 );
                ( *Bl_rel_t_minus_halfdt )( i, j )= ( *Bl_rel_t_minus_halfdt )( i, j+1 );
            }
        }
        
        else if( imode==1 ) {
            for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                ( *Bl_rel_t_plus_halfdt )( i, j ) = -( *Bl_rel_t_plus_halfdt )( i, j+1 );
                ( *Bl_rel_t_minus_halfdt )( i, j )= -( *Bl_rel_t_minus_halfdt )( i, j+1 );
            }
            
            for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
                ( *Br_rel_t_plus_halfdt )( i, j ) +=  i1*dt/dr/2.*( *ElAMrel )( i, j+1 )
                                              +			dt/dl/2.*( ( *EtAMrel )( i, j )-( *EtAMrel )( i-1, j ) );
                ( *Br_rel_t_minus_halfdt )( i, j )-=  i1*dt/dr/2.*( *ElAMrel )( i, j+1 )
                                              +			dt/dl/2.*( ( *EtAMrel )( i, j )-( *EtAMrel )( i-1, j ) );
            }
            for( unsigned int i=0; i<nl_d ; i++ ) {
                ( *Bt_rel_t_plus_halfdt )( i, j ) = -2.*i1*( *Br_rel_t_plus_halfdt )( i, j ) -( *Bt_rel_t_plus_halfdt )( i, j+1 );
                ( *Bt_rel_t_minus_halfdt )( i, j )= -2.*i1*( *Br_rel_t_minus_halfdt )( i, j )-( *Bt_rel_t_minus_halfdt )( i, j+1 );
            }
            
        } else { // modes > 1
            for( unsigned int  i=0 ; i<nl_p; i++ ) {
                ( *Bl_rel_t_plus_halfdt )( i, j ) = -( *Bl_rel_t_plus_halfdt )( i, j+1 );
                ( *Bl_rel_t_minus_halfdt )( i, j )= -( *Bl_rel_t_minus_halfdt )( i, j+1 );
            }
            for( unsigned int i=0 ; i<nl_d; i++ ) {
                ( *Br_rel_t_plus_halfdt )( i, j ) =0;
                ( *Br_rel_t_minus_halfdt )( i, j )=0;
            }
            for( unsigned int  i=0 ; i<nl_d ; i++ ) {
                ( *Bt_rel_t_plus_halfdt )( i, j ) = - ( *Bt_rel_t_plus_halfdt )( i, j+1 );
                ( *Bt_rel_t_minus_halfdt )( i, j )= - ( *Bt_rel_t_minus_halfdt )( i, j+1 );
            }
        }
    }


    // Final addition of the relativistic fields to the grid fields
    // Magnetic field Bl^(p,d)
    for( unsigned int i=0 ; i<nl_p;  i++ ) {
        for( unsigned int j=0 ; j<nr_d-1 ; j++ ) {
            // sum to the fields on grid
            ( *BlAM )( i, j )  += ( *Bl_rel_t_plus_halfdt )( i, j );
            ( *BlAM0 )( i, j ) += ( *Bl_rel_t_minus_halfdt )( i, j );
        }
    }
    
    // Magnetic field Br^(d,p)
    for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
        for( unsigned int j=0 ; j<nr_p ; j++ ) {
            // sum to the fields on grid
            ( *BrAM )( i, j )  += ( *Br_rel_t_plus_halfdt )( i, j );
            ( *BrAM0 )( i, j ) += ( *Br_rel_t_minus_halfdt )( i, j );
        }
    }
    
    // Magnetic field Bt^(d,d)
    for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
        for( unsigned int j=0 ; j<nr_d-1 ; j++ ) {
            ( *BtAM )( i, j )  += ( *Bt_rel_t_plus_halfdt )( i, j );
            ( *BtAM0 )( i, j ) += ( *Bt_rel_t_minus_halfdt )( i, j );
        }
    }
    
} // sum_rel_fields_to_em_fields

void ElectroMagnAM::sum_Poisson_fields_to_em_fields_AM( Patch *patch, Params &params, unsigned int imode )
{
    cField2D *ElAMPoisson  = static_cast<cField2D *>( El_Poisson_ );
    cField2D *ErAMPoisson  = static_cast<cField2D *>( Er_Poisson_ );
    cField2D *EtAMPoisson  = static_cast<cField2D *>( Et_Poisson_ );
    
    // E field already existing on the grid
    cField2D *ElAM  = static_cast<cField2D *>( El_[imode] );
    cField2D *ErAM  = static_cast<cField2D *>( Er_[imode] );
    cField2D *EtAM  = static_cast<cField2D *>( Et_[imode] );
  
    // El (d,p)
    for( unsigned int i=0; i<nl_d; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *ElAM )( i, j ) = ( *ElAM )( i, j ) + ( *ElAMPoisson )( i, j );
        }
    }
    
    // Er (p,d)
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_d; j++ ) {
            ( *ErAM )( i, j ) = ( *ErAM )( i, j ) + ( *ErAMPoisson )( i, j );
        }
    }
    
    // Et (p,p)
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *EtAM )( i, j ) = ( *EtAM )( i, j ) + ( *EtAMPoisson )( i, j );
        }
    }
    
} // sum_Poisson_fields_to_em_fields

void ElectroMagnAM::delete_phi_r_p_Ap( Patch *patch ){
    // delete temporary fields used for relativistic initialization
    delete phi_AM_;
    delete r_AM_;
    delete p_AM_;
    delete Ap_AM_;
}

void ElectroMagnAM::delete_relativistic_fields(Patch *patch){
    // delete temporary fields used for relativistic initialization
    delete El_rel_;
    delete Er_rel_;
    delete Et_rel_;
    delete Bl_rel_;
    delete Br_rel_;
    delete Bt_rel_;
    
    delete Bl_rel_t_plus_halfdt_;
    delete Br_rel_t_plus_halfdt_;
    delete Bt_rel_t_plus_halfdt_;
    delete Bl_rel_t_minus_halfdt_;
    delete Br_rel_t_minus_halfdt_;
    delete Bt_rel_t_minus_halfdt_;
}

void ElectroMagnAM::delete_Poisson_fields(Patch *patch){
    // delete temporary fields used for relativistic initialization
    delete El_Poisson_;
    delete Er_Poisson_;
    delete Et_Poisson_;
    
}

void ElectroMagnAM::initE( Patch *patch )
{
#ifdef _TODO_AM
#endif

    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;
    
} // initE


void ElectroMagnAM::centeringE( std::vector<double> E_Add )
{
    cField2D *El  = El_[0];
    cField2D *Er  = Er_[0];
    cField2D *Et  = Et_[0];
    
    // Centering electrostatic fields
    for( unsigned int i=0; i<nl_d; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *El )( i, j ) += E_Add[0];
        }
    }
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_d; j++ ) {
            ( *Er )( i, j ) += E_Add[1];
        }
    }
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *Et )( i, j ) += E_Add[2];
        }
    }
#ifdef _TODO_AM
#endif
    
} // centeringE

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::saveMagneticFields( bool is_spectral )
{
    if( is_spectral ) {
        for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
            Bl_m[imode]->deallocateDataAndSetTo( Bl_[imode] );
            Br_m[imode]->deallocateDataAndSetTo( Br_[imode] );
            Bt_m[imode]->deallocateDataAndSetTo( Bt_[imode] );
        }
    } else {
        for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
            // Static cast of the fields
            cField2D *Bl   = Bl_[imode];
            cField2D *Br   = Br_[imode];
            cField2D *Bt   = Bt_[imode];
            cField2D *Bl_old = Bl_m[imode];
            cField2D *Br_old = Br_m[imode];
            cField2D *Bt_old = Bt_m[imode];
            
            // Magnetic field Bl^(p,d)
            memcpy( &( ( *Bl_old )( 0, 0 ) ), &( ( *Bl )( 0, 0 ) ), nl_p*nr_d*sizeof( complex<double> ) );
            
            // Magnetic field Br^(d,p)
            memcpy( &( ( *Br_old )( 0, 0 ) ), &( ( *Br )( 0, 0 ) ), nl_d*nr_p*sizeof( complex<double> ) );
            
            // Magnetic field Bt^(d,d)
            memcpy( &( ( *Bt_old )( 0, 0 ) ), &( ( *Bt )( 0, 0 ) ), nl_d*nr_d*sizeof( complex<double> ) );
        }
    }
    
}//END saveMagneticFields


// Create a new field
Field *ElectroMagnAM::createField( string fieldname, Params& params )
{
    if( fieldname.substr( 0, 2 )=="El" ) {
        return FieldFactory::createComplex( dimPrim, 0, false, fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Er" ) {
        return FieldFactory::createComplex( dimPrim, 1, false, fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Et" ) {
        return FieldFactory::createComplex( dimPrim, 2, false, fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Bl" ) {
        return FieldFactory::createComplex( dimPrim, 0, true,  fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Br" ) {
        return FieldFactory::createComplex( dimPrim, 1, true,  fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Bt" ) {
        return FieldFactory::createComplex( dimPrim, 2, true,  fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Jl" ) {
        return FieldFactory::createComplex( dimPrim, 0, false, fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Jr" ) {
        return FieldFactory::createComplex( dimPrim, 1, false, fieldname, params );
    } else if( fieldname.substr( 0, 2 )=="Jt" ) {
        return FieldFactory::createComplex( dimPrim, 2, false, fieldname, params );
    } else if( fieldname.substr( 0, 3 )=="Rho" ) {
        return new cField2D( dimPrim, fieldname );
    } else if( fieldname.substr( 0, 9 )=="Env_A_abs" ) {
        return new Field2D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 7 )=="Env_Chi" ) {
        return new Field2D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 9 )=="Env_E_abs" ) {
        return new Field2D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 10 )=="Env_Ex_abs" ) {
        return new Field2D( dimPrim, 0, false, fieldname );
    }
    
    ERROR( "Cannot create field "<<fieldname );
    return NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::centerMagneticFields()
{
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
    
        // Static cast of the fields
        cField2D *Bl     = Bl_ [imode];
        cField2D *Br     = Br_ [imode];
        cField2D *Bt     = Bt_ [imode];
        cField2D *Bl_old = Bl_m[imode];
        cField2D *Br_old = Br_m[imode];
        cField2D *Bt_old = Bt_m[imode];
        
        // Magnetic field Bl^(p,d,d)
        for( unsigned int i=0 ; i<nl_p ; i++ ) {
            for( unsigned int j=0 ; j<nr_d ; j++ ) {
                ( *Bl_old )( i, j ) = ( ( *Bl )( i, j ) + ( *Bl_old )( i, j ) )*0.5;
            }
        }
        
        // Magnetic field Br^(d,p,d)
        for( unsigned int i=0 ; i<nl_d ; i++ ) {
            for( unsigned int j=0 ; j<nr_p ; j++ ) {
                ( *Br_old )( i, j ) = ( ( *Br )( i, j ) + ( *Br_old )( i, j ) )*0.5;
            }
        }
        
        // Magnetic field Bt^(d,d,p)
        for( unsigned int i=0 ; i<nl_d ; i++ ) {
            for( unsigned int j=0 ; j<nr_d ; j++ ) {
                ( *Bt_old )( i, j ) = ( ( *Bt )( i, j ) + ( *Bt_old )( i, j ) )*0.5;
            } // end for j
        } // end for i
        
    }
    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes)
{
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
    
        // Static cast of the fields
        cField2D *Jl     = static_cast<cField2D *>(Jl_ [imode]);
        cField2D *Jr     = static_cast<cField2D *>(Jr_ [imode]);
        cField2D *Jt     = static_cast<cField2D *>(Jt_ [imode]);

        // applying a single pass of the binomial filter along X
        if (ipass < passes[0]){
            // on Jl^(d,p) -- external points are treated by exchange. Boundary points not concerned by exchange are treated with a lower order filter.
            for( unsigned int i=0; i<nl_d-1; i++ ) {
                for( unsigned int j=isYmin*2; j<nr_p; j++ ) {
                        ( *Jl )( i, j) = ( ( *Jl )( i, j) + ( *Jl )( i+1, j) )*0.5;
                }
            }
            for( unsigned int i=nl_d-2; i>0; i-- ) {
                for( unsigned int j=isYmin*2; j<nr_p; j++ ) {
                        ( *Jl )( i, j) = ( ( *Jl )( i, j) + ( *Jl )( i-1, j) )*0.5;
                }
            }
            // Jr
            for( unsigned int i=0; i<nl_p-1; i++ ) {
                for( unsigned int j=isYmin*3; j<nr_d; j++ ) {
                        ( *Jr )( i, j) = ( ( *Jr )( i, j) + ( *Jr )( i+1, j) )*0.5;
                }
            }
            for( unsigned int i=nl_p-2; i>0; i-- ) {
                for( unsigned int j=isYmin*3; j<nr_d; j++ ) {
                        ( *Jr )( i, j) = ( ( *Jr )( i, j) + ( *Jr )( i-1, j) )*0.5;
                }
            }
            // Jt
            for( unsigned int i=0; i<nl_p-1; i++ ) {
                for( unsigned int j=isYmin*2; j<nr_p; j++ ) {
                        ( *Jt )( i, j) = ( ( *Jt )( i, j) + ( *Jt )( i+1, j) )*0.5;
                }
            }
            for( unsigned int i=nl_p-2; i>0; i-- ) {
                for( unsigned int j=isYmin*2; j<nr_p; j++ ) {
                        ( *Jt )( i, j) = ( ( *Jt )( i, j) + ( *Jt )( i-1, j) )*0.5;
                }
            }
        }

        // applying a single pass of the binomial filter along Y
        if (ipass < passes[1]){

            double qp[nr_p], qm[nr_p]; // Coefficient to retrieve proper binomial stencil along r in cylindrical geometry.
            complex<double> temp[nl_d];
            for( unsigned int j=isYmin*2+1; j<nr_p-1; j++ ) {
                qp[j] = 1.+dr*invR[j];
                qm[j] = 1.-dr*invR[j];
            }

            double qdp[nr_d], qdm[nr_d];
            for( unsigned int j=isYmin*3; j<nr_d-1; j++ ) {
                qdp[j] = 1.+dr*invRd[j];
                qdm[j] = 1.-dr*invRd[j];
            }

            //Jl
            if(isYmin){
                int j0 = 3;
                for( unsigned int i=1; i<nl_d-1; i++ ) {
                    temp[i] = ( qm[j0]*( *Jl )( i, j0-1) + qp[j0]*( *Jl )( i, j0) )*0.25; 
                }
            }
            for( unsigned int i=1; i<nl_d-1; i++ ) {
                for( unsigned int j=isYmin*2+1; j<nr_p-1; j++ ) {
                        ( *Jl )( i, j) = ( qm[j]*( *Jl )( i, j) + qp[j]*( *Jl )( i, j+1) )*0.5;
                }
            }
            for( unsigned int i=1; i<nl_d-1; i++ ) {
                for(  int j=nr_p-2; j>isYmin*2+1; j-- ) {
                        ( *Jl )( i, j) = ( ( *Jl )( i, j) + ( *Jl )( i, j-1)  )*0.5;
                }
            }
            if(isYmin){
                int j0 = 3;
                for( unsigned int i=1; i<nl_d-1; i++ ) {
                        ( *Jl )( i, j0) = 0.5*( *Jl )( i, j0) + temp[i];
                }
            }
            //Jr
            if(isYmin){
                int j0 = 4;
                for( unsigned int i=1; i<nl_p-1; i++ ) {
                    temp[i] = ( qdm[j0]*( *Jr )( i, j0-1) + qdp[j0]*( *Jr )( i, j0) )*0.25; 
                }
            }
            for( unsigned int i=1; i<nl_p-1; i++ ) {
                for( unsigned int j=isYmin*3+1; j<nr_d-1; j++ ) {
                        ( *Jr )( i, j) = ( qdp[j]*( *Jr )( i, j) + qdm[j]*( *Jr )( i, j+1) )*0.5;
                }
            }
            for( unsigned int i=1; i<nl_p-1; i++ ) {
                for(  int j=nr_d-2; j>isYmin*3+1; j-- ) {
                        ( *Jr )( i, j) = ( ( *Jr )( i, j) + ( *Jr )( i, j-1) )*0.5;
                }
            }
            if(isYmin){
                int j0 = 4;
                for( unsigned int i=1; i<nl_p-1; i++ ) {
                        ( *Jr )( i, j0) = 0.5*( *Jr )( i, j0) + temp[i];
                }
            }
            //Jt
            if(isYmin){
                int j0 = 3;
                for( unsigned int i=1; i<nl_p-1; i++ ) {
                    temp[i] = ( qm[j0]*( *Jt )( i, j0-1) + qp[j0]*( *Jt )( i, j0) )*0.25; 
                }
            }
            for( unsigned int i=1; i<nl_p-1; i++ ) {
                for( unsigned int j=isYmin*2+1; j<nr_p-1; j++ ) {
                        ( *Jt )( i, j) = ( qp[j]*( *Jt )( i, j) + qm[j]*( *Jt )( i, j+1) )*0.5;
                }
            }
            for( unsigned int i=1; i<nl_p-1; i++ ) {
                for(  int j=nr_p-2; j>isYmin*2+1; j-- ) {
                        ( *Jt )( i, j) = ( ( *Jt )( i, j) + ( *Jt )( i, j-1) )*0.5;
                }
            }
            if(isYmin){
                int j0 = 3;
                for( unsigned int i=1; i<nl_p-1; i++ ) {
                        ( *Jt )( i, j0) = 0.5*( *Jt )( i, j0) + temp[i];
                }
            }
        }

        //cField2D *tmp   = new cField2D( dimPrim, 0, false );
        //tmp->copyFrom( Jl );
        //for( unsigned int i=1; i<nl_d-1; i++ ) {
        //    for( unsigned int j=isYmin*2+1; j<nr_p-1; j++ ) {
        //        ( *Jl )( i, j ) = ( (   ( *tmp )( i+1, j-1 )+ 2.*( *tmp )( i, j-1 )+    ( *tmp )( i-1, j-1 ))*(double)(j_glob_+j-1)
        //                          + (2.*( *tmp )( i+1, j   )+ 4.*( *tmp )( i, j   )+ 2.*( *tmp )( i-1, j   ))*(double)(j_glob_+j  )
        //                          + (   ( *tmp )( i+1, j+1 )+ 2.*( *tmp )( i, j+1 )+    ( *tmp )( i-1, j+1 ))*(double)(j_glob_+j+1)
        //                          )/16.*dr*invR[j];
        //    }
        //}
        //delete tmp;

        //tmp   = new cField2D( dimPrim, 1, false );
        //tmp->copyFrom( Jr );
        //for( unsigned int i=1; i<nl_p-1; i++ ) {
        //    for( unsigned int j=isYmin*3+1; j<nr_d-1; j++ ) {
        //        ( *Jr )( i, j ) = ( (   ( *tmp )( i+1, j-1 )+ 2.*( *tmp )( i, j-1 )+    ( *tmp )( i-1, j-1 ))*(double)(j_glob_+j-1.5)
        //                          + (2.*( *tmp )( i+1, j   )+ 4.*( *tmp )( i, j   )+ 2.*( *tmp )( i-1, j   ))*(double)(j_glob_+j-0.5)
        //                          + (   ( *tmp )( i+1, j+1 )+ 2.*( *tmp )( i, j+1 )+    ( *tmp )( i-1, j+1 ))*(double)(j_glob_+j+0.5)
        //                          )/16.*dr*invRd[j];
        //    }
        //}
        //delete tmp;

        //tmp   = new cField2D( dimPrim, 2, false );
        //tmp->copyFrom( Jt );
        //for( unsigned int i=1; i<nl_p-1; i++ ) {
        //    for( unsigned int j=isYmin*2+1; j<nr_p-1; j++ ) {
        //        ( *Jt )( i, j ) = ( (   ( *tmp )( i+1, j-1 )+ 2.*( *tmp )( i, j-1 )+    ( *tmp )( i-1, j-1 ))*(double)(j_glob_+j-1)
        //                          + (2.*( *tmp )( i+1, j   )+ 4.*( *tmp )( i, j   )+ 2.*( *tmp )( i-1, j   ))*(double)(j_glob_+j  )
        //                          + (   ( *tmp )( i+1, j+1 )+ 2.*( *tmp )( i, j+1 )+    ( *tmp )( i-1, j+1 ))*(double)(j_glob_+j+1)
        //                          )/16.*dr*invR[j];
        //    }
        //}
        //delete tmp;
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::computeTotalRhoJ()
{
    for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
    
        // static cast of the total currents and densities
        cField2D *Jl     = Jl_[imode];
        cField2D *Jr     = Jr_[imode];
        cField2D *Jt     = Jt_[imode];
        cField2D *rho    = rho_AM_[imode];
        // -----------------------------------
        // Species currents and charge density
        // -----------------------------------
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        
            int ifield = imode*n_species+ispec;
            if( Jl_s[ifield] ) {
                cField2D *Jl2D_s  = Jl_s[ifield];
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    for( unsigned int j=0 ; j<nr_p ; j++ ) {
                        ( *Jl )( i, j ) += ( *Jl2D_s )( i, j );
                    }
                }
            }
            if( Jr_s[ifield] ) {
                cField2D *Jr2D_s  = Jr_s[ifield];
                for( unsigned int i=0 ; i<nl_p ; i++ )
                    for( unsigned int j=0 ; j<nr_d ; j++ ) {
                        ( *Jr )( i, j ) += ( *Jr2D_s )( i, j );
                    }
            }
            if( Jt_s[ifield] ) {
                cField2D *Jt2D_s  = Jt_s[ifield];
                for( unsigned int i=0 ; i<nl_p ; i++ )
                    for( unsigned int j=0 ; j<nr_p ; j++ ) {
                        ( *Jt )( i, j ) += ( *Jt2D_s )( i, j );
                    }
            }
            if( rho_AM_s[ifield] ) {
                cField2D *rho2D_s  = rho_AM_s[ifield];
                for( unsigned int i=0 ; i<nl_p ; i++ )
                    for( unsigned int j=0 ; j<nr_p ; j++ ) {
                        ( *rho )( i, j ) += ( *rho2D_s )( i, j );
                    }
            }
            
        }//END loop on species ispec
        
    }//END loop on mmodes
} //END computeTotalRhoJ

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total susceptibility from species susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::computeTotalEnvChi()
{ 
    // static cast of the total susceptibility
    Field2D *Env_Chi2Dcyl   = static_cast<Field2D *>( Env_Chi_);
    
    // -----------------------------------
    // Species susceptibility
    // -----------------------------------
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Env_Chi_s[ispec] ) {
            Field2D *Env_Chi2Dcyl_s  = static_cast<Field2D *>( Env_Chi_s[ispec] );
            for( unsigned int i=0 ; i<nl_p ; i++ ) {
                for( unsigned int j=0 ; j<nr_p ; j++ ) {
                    ( *Env_Chi2Dcyl )( i, j ) += ( *Env_Chi2Dcyl_s )( i, j );
                }
            }
        }
    }//END loop on species ispec


} //END computeTotalEnvChi


// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::computePoynting()
{
}

void ElectroMagnAM::applyExternalFields( Patch *patch )
{
#ifdef _TODO_AM
#endif
    int Nmodes = El_.size();
    
    Field *field;
    
    for( int imode=0; imode<Nmodes; imode++ ) {
        for( vector<ExtField>::iterator extfield=extFields.begin(); extfield!=extFields.end(); extfield++ ) {
            string name = LowerCase( extfield->field );
            if( El_[imode] && name==LowerCase( El_[imode]->name ) ) {
                field = El_[imode];
            } else if( Er_[imode] && name==LowerCase( Er_[imode]->name ) ) {
                field = Er_[imode];
            } else if( Et_[imode] && name==LowerCase( Et_[imode]->name ) ) {
                field = Et_[imode];
            } else if( Bl_[imode] && name==LowerCase( Bl_[imode]->name ) ) {
                field = Bl_[imode];
            } else if( Br_[imode] && name==LowerCase( Br_[imode]->name ) ) {
                field = Br_[imode];
            } else if( Bt_[imode] && name==LowerCase( Bt_[imode]->name ) ) {
                field = Bt_[imode];
            } else {
                field = NULL;
            }
            
            if( field ) {
                applyExternalField( field, extfield->profile, patch );
            };
        }
        Bl_m[imode]->copyFrom( Bl_[imode] );
        Br_m[imode]->copyFrom( Br_[imode] );
        Bt_m[imode]->copyFrom( Bt_[imode] );
    }
    
}

void ElectroMagnAM::applyPrescribedFields( Patch *patch, double time )
{

#ifdef _TODO_AM
#endif
    int Nmodes = El_.size();
    
    Field *field;

    for (int imode=0;imode<Nmodes;imode++){
        for( vector<ExtTimeField>::iterator extfield=extTimeFields.begin(); extfield!=extTimeFields.end(); extfield++ ) {
			string name = LowerCase( extfield->savedField->name );
			if( El_[imode] && name==LowerCase( El_[imode]->name ) ) {
				field = El_[imode];
			} else if( Er_[imode] && name==LowerCase( Er_[imode]->name ) ) {
				field = Er_[imode];
			} else if( Et_[imode] && name==LowerCase( Et_[imode]->name ) ) {
				field = Et_[imode];
			} else if( Bl_[imode] && name==LowerCase( Bl_[imode]->name ) ) {
				field = Bl_[imode];
			} else if( Br_[imode] && name==LowerCase( Br_[imode]->name ) ) {
				field = Br_[imode];
			} else if( Bt_[imode] && name==LowerCase( Bt_[imode]->name ) ) {
				field = Bt_[imode];
			} else {
				field = NULL;
			}
		
			if( field ){ 
				applyPrescribedField( field, extfield->profile, patch, time );
			}
        }
        Bl_m[imode]->copyFrom( Bl_[imode] );
        Br_m[imode]->copyFrom( Br_[imode] );
        Bt_m[imode]->copyFrom( Bt_[imode] );
    }

}

void ElectroMagnAM::applyExternalField( Field *my_field,  Profile *profile, Patch *patch )
{
     cField2D *field2D=static_cast<cField2D *>( my_field );
     
     vector<double> pos( 2 );
     pos[0]      = patch->getDomainLocalMin(0) + dl*( (field2D->isDual( 0 )?-0.5:0.) - oversize[0]);
     double pos1 = patch->getDomainLocalMin(1) + dr*( (field2D->isDual( 1 )?-0.5:0.) - oversize[1])+is_pxr*dr/2.;
     int N0 = ( int )field2D->dims()[0];
     int N1 = ( int )field2D->dims()[1];
     
     vector<Field *> xr( 2 );
     vector<unsigned int> n_space_to_create( 2 );
     n_space_to_create[0] = N0;
     n_space_to_create[1] = N1;
     
     for( unsigned int idim=0 ; idim<2 ; idim++ ) {
         xr[idim] = new Field2D( n_space_to_create );
     }
     
     for( int i=0 ; i<N0 ; i++ ) {
         pos[1] = pos1;
         for( int j=0 ; j<N1 ; j++ ) {
             for( unsigned int idim=0 ; idim<2 ; idim++ ) {
                 ( *xr[idim] )( i, j ) = pos[idim];
             }
             pos[1] += dr;
         }
         pos[0] += dl;
     }
     
     profile->complexValuesAt( xr, *field2D );
     
     for( unsigned int idim=0 ; idim<2 ; idim++ ) {
         delete xr[idim];
     }
     
     //for( auto &embc: emBoundCond ) {
     //    if( embc ) {
     //        embc->save_fields( my_field, patch );
     //    }
     //}
    
}

void ElectroMagnAM::applyPrescribedField( Field *my_field,  Profile *profile, Patch *patch, double time )
{
    cField2D *field2D=static_cast<cField2D *>( my_field );
    
    vector<double> pos( 2 );
    pos[0]      = dl*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( field2D->isDual( 0 )?-0.5:0. ) );
    double pos1 = dr*( ( double )( patch->getCellStartingGlobalIndex( 1 ) )+( field2D->isDual( 1 )?-0.5:0. ) );
    int N0 = ( int )field2D->dims()[0];
    int N1 = ( int )field2D->dims()[1];
    
    vector<Field *> xr( 2 );
    vector<unsigned int> n_space_to_create( 2 );
    n_space_to_create[0] = N0;
    n_space_to_create[1] = N1;

    for( unsigned int idim=0 ; idim<2 ; idim++ ) {
        xr[idim] = new Field2D( n_space_to_create );
    }

    for( int i=0 ; i<N0 ; i++ ) {
        pos[1] = pos1;
        for( int j=0 ; j<N1 ; j++ ) {
            for( unsigned int idim=0 ; idim<2 ; idim++ ) {
                ( *xr[idim] )( i, j ) = pos[idim];
            }
            pos[1] += dr;
        }
        pos[0] += dl;
    }

    profile->complexValuesAtTime( xr, time, *field2D );

    for( unsigned int idim=0 ; idim<2 ; idim++ ) {
        delete xr[idim];
    }

    //for( auto &embc: emBoundCond ) {
    //    if( embc ) {
    //        embc->save_fields( my_field, patch );
    //    }
    //}    
}


void ElectroMagnAM::initAntennas( Patch *patch, Params& params )
{

    // Filling the space profiles of antennas
    for( unsigned int i=0; i<antennas.size(); i++ ) {
        if( antennas[i].fieldName == "Jl" ) {
            antennas[i].field = FieldFactory::createComplex( dimPrim, 0, false, "Jl", params );
        } else if( antennas[i].fieldName == "Jr" ) {
            antennas[i].field = FieldFactory::createComplex( dimPrim, 1, false, "Jr", params );
        } else if( antennas[i].fieldName == "Jt" ) {
            antennas[i].field = FieldFactory::createComplex( dimPrim, 2, false, "Jt", params );
        } else {
            ERROR( "Antenna cannot be applied to field "<<antennas[i].fieldName );
        }
        
        if( antennas[i].field ) {
            applyExternalField( antennas[i].field, antennas[i].space_profile, patch );
        }
    }
    
}

