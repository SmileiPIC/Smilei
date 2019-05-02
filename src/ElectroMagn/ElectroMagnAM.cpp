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
            species_mode_name << vecSpecies[ispec]->name << "_mode_" << imode;
            Jl_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jl_" , vecSpecies[ispec]->name, "_mode_",  imode ).c_str(), dimPrim, params );
            Jr_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jr_" , vecSpecies[ispec]->name, "_mode_",  imode ).c_str(), dimPrim, params );
            Jt_s[imode*n_species+ispec]  = FieldFactory::createComplex( Tools::merge("Jt_" , vecSpecies[ispec]->name, "_mode_",  imode ).c_str(), dimPrim, params );
            rho_AM_s[imode*n_species+ispec] = new cField2D( ( "Rho_"+ species_mode_name.str() ).c_str(), dimPrim );
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
    }
    
    for( int ispec=0; ispec<nspecies*( int )nmodes; ispec++ ) {
        allFields.push_back( Jl_s[ispec] );
        allFields.push_back( Jr_s[ispec] );
        allFields.push_back( Jt_s[ispec] );
        allFields.push_back( rho_AM_s[ispec] );
    }
    
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
        delete Bl_m[imode];
        delete Br_m[imode];
        delete Bt_m[imode];
        
        delete Jl_[imode];
        delete Jr_[imode];
        delete Jt_[imode];
        delete rho_AM_[imode];
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
#ifdef _TODO_AM
    cField2D *rho = rho_AM_[0];
    
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
    
    phi_ = new cField2D( dimPrim );  // scalar potential
    r_   = new cField2D( dimPrim );  // residual vector
    p_   = new cField2D( dimPrim );  // direction vector
    Ap_  = new cField2D( dimPrim );  // A*p vector
    
    
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *phi_ )( i, j )   = 0.0;
            ( *r_ )( i, j )     = -( *rho )( i, j );
            ( *p_ )( i, j )     = ( *r_ )( i, j );
        }//j
    }//i
#endif
    
} // initPoisson

double ElectroMagnAM::compute_r()
{
    double rnew_dot_rnew_local( 0. );
    for( unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++ ) {
        for( unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++ ) {
            rnew_dot_rnew_local += ( *r_ )( i, j )*( *r_ )( i, j );
        }
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagnAM::compute_Ap( Patch *patch )
{
#ifdef _TODO_AM
#endif
} // compute_pAp

double ElectroMagnAM::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
#ifdef _TODO_AM
#endif
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagnAM::update_pand_r( double r_dot_r, double p_dot_Ap )
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *phi_ )( i, j ) += alpha_k * ( *p_ )( i, j );
            ( *r_ )( i, j )   -= alpha_k * ( *Ap_ )( i, j );
        }
    }
    
} // update_pand_r

void ElectroMagnAM::update_p( double rnew_dot_rnew, double r_dot_r )
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for( unsigned int i=0; i<nl_p; i++ ) {
        for( unsigned int j=0; j<nr_p; j++ ) {
            ( *p_ )( i, j ) = ( *r_ )( i, j ) + beta_k * ( *p_ )( i, j );
        }
    }
} // update_p

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
            Bl_m[imode] = Bl_[imode];
            Br_m[imode] = Br_[imode];
            Bt_m[imode] = Bt_[imode];
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
void ElectroMagnAM::binomialCurrentFilter()
{
    ERROR( "Binomial current filtering not yet implemented in AM" );
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
        //MESSAGE("c");
        // -----------------------------------
        // Species currents and charge density
        // -----------------------------------
        for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        
            int ifield = imode*n_species+ispec;
            // MESSAGE("cc");
            // MESSAGE(Jl_s.size());
            // MESSAGE(Jl_.size());
            // MESSAGE(ifield);
            if( Jl_s[ifield] ) {
                cField2D *Jl2D_s  = Jl_s[ifield];
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    //MESSAGE("here");
                    //MESSAGE(nr_p);
                    //MESSAGE(nl_p);
                    for( unsigned int j=0 ; j<nr_p ; j++ ) {
                        //MESSAGE("here i=" <<i << "  j="<<j);
                        ( *Jl )( i, j ) += ( *Jl2D_s )( i, j );
                    }
                }
            }
            //MESSAGE("or here");
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
    //MESSAGE("totalRj");
} //END computeTotalRhoJ

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total susceptibility from species susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnAM::computeTotalEnvChi()
{ } //END computeTotalEnvChi


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

    for (int imode=0;imode<Nmodes;imode++){
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
            
            if( field ){ 
                applyExternalField( field, extfield->profile, patch );
            };
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

//! Evaluating EM fields modes correctly on axis
void ElectroMagnAM::on_axis_J( bool diag_flag )
{

    if( isYmin ) {
    
        cField2D *Jl ;
        cField2D *Jr ;
        cField2D *Jt ;
        
        for( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        
            //static cast of the total currents and densities
            Jl    = Jl_[imode];
            Jr    = Jr_[imode];
            Jt    = Jt_[imode];
            
            // Set Jr below axis to zero for all modes
            for( unsigned int i=0; i<nl_p; i++ ) {
                for( unsigned int j=0; j<oversize[1]+1; j++ ) {
                    ( *Jr )( i, j ) = 0. ;
                }
            }
            // Set Jt below and on axis to zero for all modes except mode 1
            if( imode != 1 ) {
                for( unsigned int i=0; i<nl_p; i++ ) {
                    for( unsigned int j=0; j<oversize[1]+1; j++ ) {
                        ( *Jt )( i, j ) = 0. ;
                    }
                }
            } else {
                // Set Jt below axis to zero for all modes
                for( unsigned int i=0; i<nl_p; i++ ) {
                    for( unsigned int j=0; j<oversize[1]; j++ ) {
                        ( *Jt )( i, j ) = 0. ;
                    }
                    // Set Jt on axis for mode 1 for continuity equation
                    ( *Jt )( i, oversize[1] ) = - 1./3.* ( 4.* Icpx * ( *Jr )( i, oversize[1]+1 ) + ( *Jt )( i, oversize[1]+1 ) );
                }
            }
            if( imode != 0 ) {
                // Set Jx below and on axis to zero for all modes except mode 0
                for( unsigned int i=0; i<nl_d; i++ ) {
                    for( unsigned int j=0; j<oversize[1]+1; j++ ) {
                        ( *Jl )( i, j ) = 0. ;
                    }
                }
            } else {
                // Set Jx below axis to zero for all modes
                for( unsigned int i=0; i<nl_d; i++ ) {
                    for( unsigned int j=0; j<oversize[1]; j++ ) {
                        ( *Jl )( i, j ) = 0. ;
                    }
                }
                
            }
            
        }
        //if(diag_flag){
        //    for ( unsigned int imode=0 ; imode<nmodes ; imode++ ) {
        //        cField2D* rho   = rho_AM_[imode];
        //        if (imode == 0){
        //            for (unsigned int ism=0; ism < n_species; ism++){
        //                Jt    = Jt_s[ism];
        //                Jr    = Jr_s[ism];
        //                if ( ( Jt != NULL )  )
        //                    for (unsigned int i=0; i<nl_p; i++)
        //                        (*Jt)(i,oversize[1]) = 0. ;
        //                if ( ( Jr != NULL )  )
        //                    for (unsigned int i=0; i<nl_p; i++)
        //                        (*Jr)(i,oversize[1]) =  -(*Jr)(i,oversize[1]) ;
        
        //            }
        //        }
        //        else if (imode == 1){
        //            for (unsigned int i=0; i<nl_p; i++)
        //                (*rho)(i,oversize[1])= 0.;
        //            //Loop on all modes and species for J_s
        //            for (unsigned int ism=n_species; ism <  2*n_species; ism++){
        //                Jl    = Jl_s[ism];
        //                Jt    = Jt_s[ism];
        //                Jr    = Jr_s[ism];
        //                //if ( ( Jt != NULL ) && (Jr != NULL ) ) {
        //                //    for (unsigned int i=0; i<nl_p; i++)
        //                //        (*Jt)(i,oversize[1]) = - 1./3.* (4.* Icpx * (*Jr)(i,oversize[1]+1) + (*Jt)(i,oversize[1]+1));
        //                //    for (unsigned int i=0; i<nl_p; i++)
        //                //        (*Jr)(i,oversize[1])= 2.*Icpx* (*Jt)(i,oversize[1])-(*Jr)(i,oversize[1]+1) ;
        //                //}
        //                if ( Jl != NULL )
        //                    for (unsigned int i=0; i<nl_d; i++)
        //                        (*Jl)(i,oversize[1])= 0. ;
        //            }
        //        }
        //        else {  // imode > 1
        //            //Loop on all modes and species for J_s
        //            for (unsigned int ism=2*n_species; ism <  n_species*nmodes ; ism++){
        //                Jt    = Jt_s[ism];
        //                if ( Jt != NULL )
        //                    for (unsigned int i=0; i<nl_p; i++)
        //                        (*Jt)(i,oversize[1]) = 0.;
        
        //                Jr    = Jr_s[ism];
        //                if ( Jr != NULL )
        //                    for (unsigned int i=0; i<nl_p; i++)
        //                        (*Jr)(i,oversize[1]) = -(*Jr)(i,oversize[1]+1) ;
        
        //                Jl    = Jl_s[ism];
        //                if ( Jl != NULL )
        //                    for (unsigned int i=0; i<nl_d; i++)
        //                        (*Jl)(i,oversize[1])= 0. ;
        //            }
        //        }
        //    }
        //}
    }
    return;
}

