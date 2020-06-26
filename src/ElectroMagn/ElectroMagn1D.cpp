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
ElectroMagn1D::ElectroMagn1D( Params &params, DomainDecomposition *domain_decomposition, vector<Species *> &vecSpecies, Patch *patch )
    : ElectroMagn( params, domain_decomposition, vecSpecies, patch )
{
    initElectroMagn1DQuantities( params, patch );
    
    // Charge and current densities for each species
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        Jx_s[ispec]  = new Field1D( Tools::merge( "Jx_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        Jy_s[ispec]  = new Field1D( Tools::merge( "Jy_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        Jz_s[ispec]  = new Field1D( Tools::merge( "Jz_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        rho_s[ispec] = new Field1D( Tools::merge( "Rho_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        
        if( params.Laser_Envelope_model ) {
            Env_Chi_s[ispec] = new Field1D( Tools::merge( "Env_Chi_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        }
    }
    
}//END constructor Electromagn1D


ElectroMagn1D::ElectroMagn1D( ElectroMagn1D *emFields, Params &params, Patch *patch )
    : ElectroMagn( emFields, params, patch )
{
    initElectroMagn1DQuantities( params, patch );
    
    // Charge and current densities for each species
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( emFields->Jx_s[ispec] != NULL ) {
            if( emFields->Jx_s[ispec]->data_ != NULL ) {
                Jx_s[ispec]  = new Field1D( dimPrim, 0, false, emFields->Jx_s[ispec]->name );
            } else {
                Jx_s[ispec]  = new Field1D( emFields->Jx_s[ispec]->name, dimPrim );
            }
        }
        if( emFields->Jy_s[ispec] != NULL ) {
            if( emFields->Jy_s[ispec]->data_ != NULL ) {
                Jy_s[ispec]  = new Field1D( dimPrim, 1, false, emFields->Jy_s[ispec]->name );
            } else {
                Jy_s[ispec]  = new Field1D( emFields->Jy_s[ispec]->name, dimPrim );
            }
        }
        if( emFields->Jz_s[ispec] != NULL ) {
            if( emFields->Jz_s[ispec]->data_ != NULL ) {
                Jz_s[ispec]  = new Field1D( dimPrim, 2, false, emFields->Jz_s[ispec]->name );
            } else {
                Jz_s[ispec]  = new Field1D( emFields->Jz_s[ispec]->name, dimPrim );
            }
        }
        if( emFields->rho_s[ispec] != NULL ) {
            if( emFields->rho_s[ispec]->data_ != NULL ) {
                rho_s[ispec] = new Field1D( dimPrim, emFields->rho_s[ispec]->name );
            } else {
                rho_s[ispec]  = new Field1D( emFields->rho_s[ispec]->name, dimPrim );
            }
        }
        
        if( params.Laser_Envelope_model ) {
            if( emFields->Env_Chi_s[ispec] != NULL ) {
                if( emFields->Env_Chi_s[ispec]->data_ != NULL ) {
                    Env_Chi_s[ispec] = new Field1D( dimPrim, emFields->Env_Chi_s[ispec]->name );
                } else {
                    Env_Chi_s[ispec]  = new Field1D( emFields->Env_Chi_s[ispec]->name, dimPrim );
                }
            }
        }
        
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn1D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::initElectroMagn1DQuantities( Params &params, Patch *patch )
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
    for( size_t i=0 ; i<nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2;
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    
    // Allocation of the EM fields
    Ex_  = new Field1D( dimPrim, 0, false, "Ex" );
    Ey_  = new Field1D( dimPrim, 1, false, "Ey" );
    Ez_  = new Field1D( dimPrim, 2, false, "Ez" );
    Bx_  = new Field1D( dimPrim, 0, true,  "Bx" );
    By_  = new Field1D( dimPrim, 1, true,  "By" );
    Bz_  = new Field1D( dimPrim, 2, true,  "Bz" );
    Bx_m = new Field1D( dimPrim, 0, true,  "Bx_m" );
    By_m = new Field1D( dimPrim, 1, true,  "By_m" );
    Bz_m = new Field1D( dimPrim, 2, true,  "Bz_m" );
    
    if( params.Laser_Envelope_model ) {
        Env_A_abs_  = new Field1D( dimPrim, "Env_A_abs" );
        Env_Chi_    = new Field1D( dimPrim, "Env_Chi" );
        Env_E_abs_  = new Field1D( dimPrim, "Env_E_abs" );
        Env_Ex_abs_ = new Field1D( dimPrim, "Env_Ex_abs" );
    }
    // Total charge currents and densities
    Jx_   = new Field1D( dimPrim, 0, false, "Jx" );
    Jy_   = new Field1D( dimPrim, 1, false, "Jy" );
    Jz_   = new Field1D( dimPrim, 2, false, "Jz" );
    rho_  = new Field1D( dimPrim, "Rho" );
    
    
    // ----------------------------------------------------------------
    // Definition of the min and max index according to chosen oversize
    // ----------------------------------------------------------------
    index_bc_min.resize( nDim_field, 0 );
    index_bc_max.resize( nDim_field, 0 );
    for( size_t i=0 ; i<nDim_field ; i++ ) {
        index_bc_min[i] = oversize[i];
        index_bc_max[i] = dimDual[i]-oversize[i]-1;
    }
    
    
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
                    if( ( patch->Pcoordinates[i]!=0 ) && ( patch->Pcoordinates[i]!=( unsigned int )params.number_of_patches[i]-1 ) ) {
                        bufsize[i][isDual]--;
                    }
                }
                
            } // if ( params.number_of_patches[i]!=1 )
        } // for (int isDual=0 ; isDual
    } // for (unsigned int i=0 ; i<nDim_field
}

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

void ElectroMagn1D::initPoisson( Patch *patch )
{
    Field1D *rho1D = static_cast<Field1D *>( rho_ );
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------
    
    index_min_p_.resize( 1, 0 );
    index_max_p_.resize( 1, 0 );
    
    index_min_p_[0] = oversize[0];
    index_max_p_[0] = nx_p - 2 - oversize[0];
    if( patch->isXmin() ) {
        index_min_p_[0] = 0;
    }
    if( patch->isXmax() ) {
        index_max_p_[0] = nx_p-1;
    }
    
    phi_ = new Field1D( dimPrim );  // scalar potential
    r_   = new Field1D( dimPrim );  // residual vector
    p_   = new Field1D( dimPrim );  // direction vector
    Ap_  = new Field1D( dimPrim );  // A*p vector
    
    // double       dx_sq          = dx*dx;
    
    // phi: scalar potential, r: residual and p: direction
    for( unsigned int i=0 ; i<dimPrim[0] ; i++ ) {
        ( *phi_ )( i )   = 0.0;
        //(*r_)(i)     = -dx_sq * (*rho1D)(i);
        ( *r_ )( i )     = - ( *rho1D )( i );
        ( *p_ )( i )     = ( *r_ )( i );
    }
} // initPoisson

double ElectroMagn1D::compute_r()
{
    double rnew_dot_rnew_local( 0. );
    for( unsigned int i=index_min_p_[0] ; i<=index_max_p_[0] ; i++ ) {
        rnew_dot_rnew_local += ( *r_ )( i )*( *r_ )( i );
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn1D::compute_Ap( Patch *patch )
{

    double one_ov_dx_sq       = 1.0/( dx*dx );
    double two_ov_dx2         = 2.0*( 1.0/( dx*dx ) );
    
    // vector product Ap = A*p
    for( unsigned int i=1 ; i<dimPrim[0]-1 ; i++ ) {
        ( *Ap_ )( i ) = one_ov_dx_sq * ( ( *p_ )( i-1 ) + ( *p_ )( i+1 ) )  - two_ov_dx2*( *p_ )( i )   ;
    }
    
    // apply BC on Ap
    if( patch->isXmin() ) {
        ( *Ap_ )( 0 )      = one_ov_dx_sq * ( ( *p_ )( 1 ) )      - two_ov_dx2*( *p_ )( 0 );
    }
    if( patch->isXmax() ) {
        ( *Ap_ )( nx_p-1 ) = one_ov_dx_sq * ( ( *p_ )( nx_p-2 ) ) - two_ov_dx2*( *p_ )( nx_p-1 );
    }
    
} // compute_Ap

void ElectroMagn1D::compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean )
{

    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    double one_ov_dx_sq_ov_gamma_sq       = 1.0/( dx*dx )/( gamma_mean*gamma_mean );
    double two_ov_dxgam2                  = 2.0*( 1.0/( dx*dx )/( gamma_mean*gamma_mean ) );
    
    // vector product Ap = A*p
    for( unsigned int i=1 ; i<dimPrim[0]-1 ; i++ ) {
        ( *Ap_ )( i ) = one_ov_dx_sq_ov_gamma_sq * ( ( *p_ )( i-1 ) + ( *p_ )( i+1 ) ) - two_ov_dxgam2 *( *p_ )( i )   ;
    }
    
    // apply BC on Ap
    if( patch->isXmin() ) {
        ( *Ap_ )( 0 )      = one_ov_dx_sq_ov_gamma_sq * ( ( *p_ )( 1 ) )     - two_ov_dxgam2 * ( *p_ )( 0 );
    }
    if( patch->isXmax() ) {
        ( *Ap_ )( nx_p-1 ) = one_ov_dx_sq_ov_gamma_sq * ( ( *p_ )( nx_p-2 ) )- two_ov_dxgam2 * ( *p_ )( nx_p-1 );
    }
    
} // compute_Ap_relativistic_Poisson

double ElectroMagn1D::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    for( unsigned int i=index_min_p_[0] ; i<=index_max_p_[0] ; i++ ) {
        p_dot_Ap_local += ( *p_ )( i )*( *Ap_ )( i );
    }
    return p_dot_Ap_local;
    
} // compute_pAp

void ElectroMagn1D::update_pand_r( double r_dot_r, double p_dot_Ap )
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for( unsigned int i=0 ; i<dimPrim[0] ; i++ ) {
        ( *phi_ )( i ) += alpha_k * ( *p_ )( i );
        ( *r_ )( i )   -= alpha_k * ( *Ap_ )( i );
    }
    
} // update_pand_r

void ElectroMagn1D::update_p( double rnew_dot_rnew, double r_dot_r )
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for( unsigned int i=0 ; i<dimPrim[0] ; i++ ) {
        ( *p_ )( i ) = ( *r_ )( i ) + beta_k * ( *p_ )( i );
    }
} // update_p

void ElectroMagn1D::initE( Patch *patch )
{
    Field1D *Ex1D  = static_cast<Field1D *>( Ex_ );
    Field1D *rho1D = static_cast<Field1D *>( rho_ );
    
    // ----------------------------------
    // Compute the electrostatic field Ex
    // ----------------------------------
    
    for( unsigned int i=1; i<nx_d-1; i++ ) {
        ( *Ex1D )( i ) = ( ( *phi_ )( i-1 )-( *phi_ )( i ) )/dx;
    }
    
    // BC on Ex
    if( patch->isXmin() ) {
        ( *Ex1D )( 0 )      = ( *Ex1D )( 1 )      - dx*( *rho1D )( 0 );
    }
    if( patch->isXmax() ) {
        ( *Ex1D )( nx_d-1 ) = ( *Ex1D )( nx_d-2 ) + dx*( *rho1D )( nx_p-1 );
    }
    
    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;
    
} // initE

void ElectroMagn1D::initE_relativistic_Poisson( Patch *patch, double gamma_mean )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    Field1D *Ex1D  = static_cast<Field1D *>( Ex_rel_ );
    Field1D *rho1D = static_cast<Field1D *>( rho_ );
    
    // ----------------------------------
    // Compute the electrostatic field Ex
    // ----------------------------------
    
    
    
    for( unsigned int i=1; i<nx_p-1; i++ ) {
        ( *Ex1D )( i ) = ( ( *phi_ )( i-1 )-( *phi_ )( i ) )/dx/gamma_mean/gamma_mean;
    }
    
    // BC on Ex
    if( patch->isXmin() ) {
        ( *Ex1D )( 0 )      = ( *Ex1D )( 1 )      - dx*( *rho1D )( 0 );
    }
    if( patch->isXmax() ) {
        ( *Ex1D )( nx_d-1 ) = ( *Ex1D )( nx_d-2 ) + dx*( *rho1D )( nx_p-1 );
    }
    
    
    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;
    
} // initE_relativistic_Poisson

void ElectroMagn1D::initB_relativistic_Poisson( Patch *patch, double gamma_mean )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    // For some inconsistency the B field in 1D seems zero - am I wrong?
    
} // initB_relativistic_Poisson

void ElectroMagn1D::center_fields_from_relativistic_Poisson( Patch *patch )
{

    // In 1D no centering is necessary, as E is already centered and there is no field B in relativistic initialization
    
}

void ElectroMagn1D::initRelativisticPoissonFields( Patch *patch )
{
    // init temporary fields for relativistic field initialization, to be added to the already present electromagnetic fields
    
    Ex_rel_  = new Field1D( dimPrim, 0, false, "Ex_rel" );
    Ey_rel_  = new Field1D( dimPrim, 1, false, "Ey_rel" );
    Ez_rel_  = new Field1D( dimPrim, 2, false, "Ez_rel" );
    Bx_rel_  = new Field1D( dimPrim, 0, true,  "Bx_rel" ); // will be identically zero
    By_rel_  = new Field1D( dimPrim, 1, false,  "By_rel" ); // is equal to -beta*Ez, thus it inherits the same centering of Ez
    Bz_rel_  = new Field1D( dimPrim, 2, false,  "Bz_rel" ); // is equal to  beta*Ey, thus it inherits the same centering of Ey
    
} // initRelativisticPoissonFields

void ElectroMagn1D::sum_rel_fields_to_em_fields( Patch *patch )
{
    Field1D *Ex1Drel  = static_cast<Field1D *>( Ex_rel_ );
    Field1D *Ex1D  = static_cast<Field1D *>( Ex_ );
    
    
    // Ex
    for( unsigned int i=0; i<nx_d; i++ ) {
        ( *Ex1D )( i ) = ( *Ex1D )( i ) + ( *Ex1Drel )( i );
    }
    
    // delete temporary fields used for relativistic initialization
    delete Ex_rel_;
    delete Ey_rel_;
    delete Ez_rel_;
    delete Bx_rel_;
    delete By_rel_;
    delete Bz_rel_;
    
    
} // sum_rel_fields_to_em_fields


void ElectroMagn1D::centeringE( std::vector<double> E_Add )
{
    Field1D *Ex1D  = static_cast<Field1D *>( Ex_ );
    for( unsigned int i=0; i<nx_d; i++ ) {
        ( *Ex1D )( i ) += E_Add[0];
    }
    
} // centeringE

void ElectroMagn1D::centeringErel( std::vector<double> E_Add )
{
    Field1D *Ex1D  = static_cast<Field1D *>( Ex_rel_ );
    for( unsigned int i=0; i<nx_d; i++ ) {
        ( *Ex1D )( i ) += E_Add[0];
    }
    
} // centeringErel

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::saveMagneticFields( bool is_spectral )
{
    // Static cast of the fields
    if( !is_spectral ) {
        Field1D *Bx1D   = static_cast<Field1D *>( Bx_ );
        Field1D *By1D   = static_cast<Field1D *>( By_ );
        Field1D *Bz1D   = static_cast<Field1D *>( Bz_ );
        Field1D *Bx1D_m = static_cast<Field1D *>( Bx_m );
        Field1D *By1D_m = static_cast<Field1D *>( By_m );
        Field1D *Bz1D_m = static_cast<Field1D *>( Bz_m );
        
        // for Bx^(p)
        for( unsigned int i=0 ; i<dimPrim[0] ; i++ ) {
            ( *Bx1D_m )( i )=( *Bx1D )( i );
        }
        //for By^(d) & Bz^(d)
        for( unsigned int i=0 ; i<dimDual[0] ; i++ ) {
            ( *By1D_m )( i ) = ( *By1D )( i );
            ( *Bz1D_m )( i ) = ( *Bz1D )( i );
        }
    } else {
        Bx_m->deallocateDataAndSetTo( Bx_ );
        By_m->deallocateDataAndSetTo( By_ );
        Bz_m->deallocateDataAndSetTo( Bz_ );
    }
    
}//END saveMagneticFields



//// ---------------------------------------------------------------------------------------------------------------------
//// Maxwell solver using the FDTD scheme
//// ---------------------------------------------------------------------------------------------------------------------
//void ElectroMagn1D::solveMaxwellAmpere()
//{
//
//    Field1D* Ex1D = static_cast<Field1D*>(Ex_);
//    Field1D* Ey1D = static_cast<Field1D*>(Ey_);
//    Field1D* Ez1D = static_cast<Field1D*>(Ez_);
//    Field1D* By1D = static_cast<Field1D*>(By_);
//    Field1D* Bz1D = static_cast<Field1D*>(Bz_);
//    Field1D* Jx1D = static_cast<Field1D*>(Jx_);
//    Field1D* Jy1D = static_cast<Field1D*>(Jy_);
//    Field1D* Jz1D = static_cast<Field1D*>(Jz_);
//
//    // --------------------
//    // Solve Maxwell-Ampere
//    // --------------------
//    // Calculate the electrostatic field ex on the dual grid
//    //for (unsigned int ix=0 ; ix<nx_d ; ix++){
//    for (unsigned int ix=0 ; ix<dimDual[0] ; ix++) {
//        (*Ex1D)(ix)= (*Ex1D)(ix) - timestep* (*Jx1D)(ix) ;
//    }
//    // Transverse fields ey, ez  are defined on the primal grid
//    //for (unsigned int ix=0 ; ix<nx_p ; ix++) {
//    for (unsigned int ix=0 ; ix<dimPrim[0] ; ix++) {
//        (*Ey1D)(ix)= (*Ey1D)(ix) - dt_ov_dx * ( (*Bz1D)(ix+1) - (*Bz1D)(ix)) - timestep * (*Jy1D)(ix) ;
//        (*Ez1D)(ix)= (*Ez1D)(ix) + dt_ov_dx * ( (*By1D)(ix+1) - (*By1D)(ix)) - timestep * (*Jz1D)(ix) ;
//    }
//
//}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::centerMagneticFields()
{
    // Static cast of the fields
    Field1D *Bx1D   = static_cast<Field1D *>( Bx_ );
    Field1D *By1D   = static_cast<Field1D *>( By_ );
    Field1D *Bz1D   = static_cast<Field1D *>( Bz_ );
    Field1D *Bx1D_m = static_cast<Field1D *>( Bx_m );
    Field1D *By1D_m = static_cast<Field1D *>( By_m );
    Field1D *Bz1D_m = static_cast<Field1D *>( Bz_m );
    
    // for Bx^(p)
    for( unsigned int i=0 ; i<dimPrim[0] ; i++ ) {
        ( *Bx1D_m )( i ) = ( ( *Bx1D )( i )+ ( *Bx1D_m )( i ) )*0.5 ;
    }
    
    // for By^(d) & Bz^(d)
    for( unsigned int i=0 ; i<dimDual[0] ; i++ ) {
        ( *By1D_m )( i )= ( ( *By1D )( i )+( *By1D_m )( i ) )*0.5 ;
        ( *Bz1D_m )( i )= ( ( *Bz1D )( i )+( *Bz1D_m )( i ) )*0.5 ;
    }
    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes)
{

    Field1D *Jx1D     = static_cast<Field1D *>( Jx_ );
    Field1D *Jy1D     = static_cast<Field1D *>( Jy_ );
    Field1D *Jz1D     = static_cast<Field1D *>( Jz_ );
    
    // Apply a single pass of the binomial filter on currents
    
    // on Jx^(d) -- external points are treated by exchange
    double temp0, tempN;
    //Field1D *tmp   = new Field1D(dimPrim, 0, false);
    //tmp->copyFrom(Jx1D);
    temp0 = ( *Jx1D )( 0 );
    tempN = ( *Jx1D )( dimDual[0]-1 );
    for( unsigned int ix=0 ; ix<dimDual[0]-1 ; ix++ ) {
        ( *Jx1D )( ix )  = ( ( *Jx1D )( ix ) + ( *Jx1D )( ix+1 ) ) * 0.5 ;
    }
    for( unsigned int ix=1 ; ix<dimDual[0]-1 ; ix++ ) {
        ( *Jx1D )( ix )  = ( ( *Jx1D )( ix-1 ) + ( *Jx1D )( ix ) ) * 0.5 ;
    }
    ( *Jx1D )( 0 ) = temp0;
    ( *Jx1D )( dimDual[0]-1 ) = tempN;
    
    temp0 = ( *Jy1D )( 0 );
    tempN = ( *Jy1D )( dimPrim[0]-1 );
    // on Jy^(p) -- external points are treated by exchange
    //tmp   = new Field1D(dimPrim, 1, false);
    //tmp->copyFrom(Jy1D);
    for( unsigned int ix=0 ; ix<dimPrim[0]-1 ; ix++ ) {
        ( *Jy1D )( ix )  = ( ( *Jy1D )( ix ) + ( *Jy1D )( ix+1 ) ) * 0.5 ;
    }
    for( unsigned int ix=1 ; ix<dimPrim[0]-1 ; ix++ ) {
        ( *Jy1D )( ix )  = ( ( *Jy1D )( ix ) + ( *Jy1D )( ix+1 ) ) * 0.5 ;
    }
    //delete tmp;
    
    ( *Jy1D )( 0 ) = temp0;
    ( *Jy1D )( dimPrim[0]-1 ) = tempN;
    
    temp0 = ( *Jz1D )( 0 );
    tempN = ( *Jz1D )( dimPrim[0]-1 );
    // on Jz^(p) -- external points are treated by exchange
    //tmp   = new Field1D(dimPrim, 2, false);
    //tmp->copyFrom(Jz1D);
    for( unsigned int ix=0 ; ix<dimPrim[0]-1 ; ix++ ) {
        ( *Jz1D )( ix )  = ( ( *Jz1D )( ix ) + ( *Jz1D )( ix+1 ) ) * 0.5 ;
    }
    for( unsigned int ix=1 ; ix<dimPrim[0]-1 ; ix++ ) {
        ( *Jz1D )( ix )  = ( ( *Jz1D )( ix ) + ( *Jz1D )( ix+1 ) ) * 0.5 ;
    }
    //delete tmp;
    ( *Jz1D )( 0 ) = temp0;
    ( *Jz1D )( dimPrim[0]-1 ) = tempN;
    
}

void ElectroMagn1D::customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes, std::vector<double> filtering_coeff)
{
    // Static-cast of the currents
    Field1D *Jx1D = static_cast<Field1D *>( Jx_ );
    Field1D *Jy1D = static_cast<Field1D *>( Jy_ );
    Field1D *Jz1D = static_cast<Field1D *>( Jz_ );

    // Upsampling factor
    // m=1 : No upsampling ......................... f_Nyquist *= 1
    // m=2 : 1 zero(s) between two data points ..... f_Nyquist *= 2
    // m=3 : 2 zero(s) between two data points ..... f_Nyquist *= 3
    // m=4 : 3 zero(s) between two data points ..... f_Nyquist *= 4
    unsigned int m=1 ;

    // Guard-Cell Current
    unsigned int gcfilt=0 ;

    // Applying a single pass of the custom FIR based filter along X
    if (ipass < passes[0]){
        Field1D *tmp   = new Field1D( dimPrim, 0, false );
        tmp->copyFrom( Jx1D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_d-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            ( *Jx1D )( i ) = 0. ;
            for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                ( *Jx1D )( i ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
            }
            ( *Jx1D )( i ) *= m ;
        }
        delete tmp;
        tmp   = new Field1D( dimPrim, 1, false );
        tmp->copyFrom( Jy1D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            ( *Jy1D )( i ) = 0. ;
            for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                ( *Jy1D )( i ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
            }
            ( *Jy1D )( i ) *= m ;
        }
        delete tmp;
        tmp   = new Field1D( dimPrim, 2, false );
        tmp->copyFrom( Jz1D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            ( *Jz1D )( i ) = 0. ;
            for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
               ( *Jz1D )( i ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
            }
            ( *Jz1D )( i ) *= m ;
        }
        delete tmp;
    }
}



// Create a new field
Field *ElectroMagn1D::createField( string fieldname, Params& params )
{
    if( fieldname.substr( 0, 2 )=="Ex" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 2 )=="Ey" ) {
        return new Field1D( dimPrim, 1, false, fieldname );
    } else if( fieldname.substr( 0, 2 )=="Ez" ) {
        return new Field1D( dimPrim, 2, false, fieldname );
    } else if( fieldname.substr( 0, 2 )=="Bx" ) {
        return new Field1D( dimPrim, 0, true,  fieldname );
    } else if( fieldname.substr( 0, 2 )=="By" ) {
        return new Field1D( dimPrim, 1, true,  fieldname );
    } else if( fieldname.substr( 0, 2 )=="Bz" ) {
        return new Field1D( dimPrim, 2, true,  fieldname );
    } else if( fieldname.substr( 0, 2 )=="Jx" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 2 )=="Jy" ) {
        return new Field1D( dimPrim, 1, false, fieldname );
    } else if( fieldname.substr( 0, 2 )=="Jz" ) {
        return new Field1D( dimPrim, 2, false, fieldname );
    } else if( fieldname.substr( 0, 3 )=="Rho" ) {
        return new Field1D( dimPrim, fieldname );
    } else if( fieldname.substr( 0, 9 )=="Env_A_abs" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 7 )=="Env_Chi" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 9 )=="Env_E_abs" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    } else if( fieldname.substr( 0, 10 )=="Env_Ex_abs" ) {
        return new Field1D( dimPrim, 0, false, fieldname );
    }
    
    ERROR( "Cannot create field "<<fieldname );
    return NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalRhoJ()
{
    Field1D *Jx1D    = static_cast<Field1D *>( Jx_ );
    Field1D *Jy1D    = static_cast<Field1D *>( Jy_ );
    Field1D *Jz1D    = static_cast<Field1D *>( Jz_ );
    Field1D *rho1D   = static_cast<Field1D *>( rho_ );
    
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Jx_s[ispec] ) {
            Field1D *Jx1D_s  = static_cast<Field1D *>( Jx_s[ispec] );
            for( unsigned int ix=0 ; ix<=dimPrim[0] ; ix++ ) {
                ( *Jx1D )( ix )  += ( *Jx1D_s )( ix );
            }
        }
        if( Jy_s[ispec] ) {
            Field1D *Jy1D_s  = static_cast<Field1D *>( Jy_s[ispec] );
            for( unsigned int ix=0 ; ix<dimPrim[0] ; ix++ ) {
                ( *Jy1D )( ix )  += ( *Jy1D_s )( ix );
            }
        }
        if( Jz_s[ispec] ) {
            Field1D *Jz1D_s  = static_cast<Field1D *>( Jz_s[ispec] );
            for( unsigned int ix=0 ; ix<dimPrim[0] ; ix++ ) {
                ( *Jz1D )( ix )  += ( *Jz1D_s )( ix );
            }
        }
        if( rho_s[ispec] ) {
            Field1D *rho1D_s  = static_cast<Field1D *>( rho_s[ispec] );
            for( unsigned int ix=0 ; ix<dimPrim[0] ; ix++ ) {
                ( *rho1D )( ix )  += ( *rho1D_s )( ix );
            }
        }
    }//END loop on species ispec
}

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total susceptibility from species susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn1D::computeTotalEnvChi()
{
    // static cast of the total susceptibility
    Field1D *Env_Chi1D   = static_cast<Field1D *>( Env_Chi_ );
    
    // -----------------------------------
    // Species susceptibility
    // -----------------------------------
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Env_Chi_s[ispec] ) {
            Field1D *Env_Chi1D_s  = static_cast<Field1D *>( Env_Chi_s[ispec] );
            for( unsigned int i=0 ; i<nx_p ; i++ ) {
                ( *Env_Chi1D )( i ) += ( *Env_Chi1D_s )( i );
            }
        }
    }//END loop on species ispec
    
} // END computeTotalEnvChi
// --------------------------------------------------------------------------
// Compute Poynting (return the electromagnetic energy injected at the border
// --------------------------------------------------------------------------
void ElectroMagn1D::computePoynting()
{

    // Xmin border (Energy injected = +Poynting)
    if( isXmin ) {
        unsigned int iEy=istart[0][Ey_->isDual( 0 )];
        unsigned int iBz=istart[0][Bz_m->isDual( 0 )];
        unsigned int iEz=istart[0][Ez_->isDual( 0 )];
        unsigned int iBy=istart[0][By_m->isDual( 0 )];
        
        poynting_inst[0][0]=0.5*timestep*( ( *Ey_ )( iEy ) * ( ( *Bz_m )( iBz ) + ( *Bz_m )( iBz+1 ) ) -
                                           ( *Ez_ )( iEz ) * ( ( *By_m )( iBy ) + ( *By_m )( iBy+1 ) ) );
        poynting[0][0] += poynting_inst[0][0];
    }
    
    // Xmax border (Energy injected = -Poynting)
    if( isXmax ) {
        unsigned int offset = bufsize[0][Ey_->isDual( 0 )];
        
        unsigned int iEy=istart[0][Ey_ ->isDual( 0 )] + offset;
        unsigned int iBz=istart[0][Bz_m->isDual( 0 )] + offset;
        unsigned int iEz=istart[0][Ez_ ->isDual( 0 )] + offset;
        unsigned int iBy=istart[0][By_m->isDual( 0 )] + offset;
        
        poynting_inst[1][0]=0.5*timestep*( ( *Ey_ )( iEy ) * ( ( *Bz_m )( iBz ) + ( *Bz_m )( iBz+1 ) ) -
                                           ( *Ez_ )( iEz ) * ( ( *By_m )( iBy ) + ( *By_m )( iBy+1 ) ) );
        poynting[1][0] -= poynting_inst[1][0];
        
    }
}

void ElectroMagn1D::applyExternalField( Field *my_field,  Profile *profile, Patch *patch )
{
    Field1D *field1D=static_cast<Field1D *>( my_field );
    if( patch->hindex==0 ) {
        MESSAGE( my_field->name );
    }
    
    vector<double> pos( 1 );
    pos[0] = dx * ( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( field1D->isDual( 0 )?-0.5:0. ) );
    int N = ( int )field1D->dims()[0];
    
    // USING UNSIGNED INT CREATES PB WITH PERIODIC BCs
    for( int i=0 ; i<N ; i++ ) {
        ( *field1D )( i ) += profile->valueAt( pos );
        pos[0] += dx;
    }
    
}

void ElectroMagn1D::applyPrescribedField( Field *my_field,  Profile *profile, Patch *patch, double time )
{
    Field1D *field1D=static_cast<Field1D *>( my_field );
    
    vector<double> pos( 1 );
    pos[0] = dx * ( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( field1D->isDual( 0 )?-0.5:0. ) );
    int N = ( int )field1D->dims()[0];
    
    // USING UNSIGNED INT CREATES PB WITH PERIODIC BCs
    for( int i=0 ; i<N ; i++ ) {
        ( *field1D )( i ) += profile->valueAt( pos, time ); 
        pos[0] += dx;
    }
    
}



void ElectroMagn1D::initAntennas( Patch *patch, Params& params )
{

    // Filling the space profiles of antennas
    for( unsigned int i=0; i<antennas.size(); i++ ) {
        if( antennas[i].fieldName == "Jx" ) {
            antennas[i].field = new Field1D( dimPrim, 0, false, "Jx" );
        } else if( antennas[i].fieldName == "Jy" ) {
            antennas[i].field = new Field1D( dimPrim, 1, false, "Jy" );
        } else if( antennas[i].fieldName == "Jz" ) {
            antennas[i].field = new Field1D( dimPrim, 2, false, "Jz" );
        }
        
        if( antennas[i].field ) {
            applyExternalField( antennas[i].field, antennas[i].space_profile, patch );
        }
    }
    
}
