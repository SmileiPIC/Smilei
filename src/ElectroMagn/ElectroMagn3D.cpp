#include "ElectroMagn3D.h"

#include <cmath>

#include <iostream>
#include <sstream>

#include "Params.h"
#include "Field3D.h"
#include "FieldFactory.h"

#include "Patch.h"
#include <cstring>

#include "Profile.h"

#include "ElectroMagnBC.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Electromagn3D
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagn3D::ElectroMagn3D( Params &params, DomainDecomposition *domain_decomposition, vector<Species *> &vecSpecies, Patch *patch ) :
    ElectroMagn( params, domain_decomposition, vecSpecies, patch ),
    isYmin( patch->isYmin() ),
    isYmax( patch->isYmax() ),
    isZmax( patch->isZmax() ),
    isZmin( patch->isZmin() )
{

    initElectroMagn3DQuantities( params, patch );
    
    // Charge currents currents and density for each species
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        Jx_s[ispec]  = FieldFactory::create( Tools::merge("Jx_" , vecSpecies[ispec]->name_ ).c_str(), dimPrim, params );
        Jy_s[ispec]  = FieldFactory::create( Tools::merge("Jy_" , vecSpecies[ispec]->name_ ).c_str(), dimPrim, params );
        Jz_s[ispec]  = FieldFactory::create( Tools::merge("Jz_" , vecSpecies[ispec]->name_ ).c_str(), dimPrim, params );
        rho_s[ispec] = new Field3D( Tools::merge( "Rho_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        
        if( params.Laser_Envelope_model ) {
            Env_Chi_s[ispec] = new Field3D( Tools::merge( "Env_Chi_", vecSpecies[ispec]->name_ ).c_str(), dimPrim );
        }
        
    }
    
}//END constructor Electromagn3D


ElectroMagn3D::ElectroMagn3D( ElectroMagn3D *emFields, Params &params, Patch *patch ) :
    ElectroMagn( emFields, params, patch ),
    isYmin( patch->isYmin() ),
    isYmax( patch->isYmax() ),
    isZmax( patch->isZmax() ),
    isZmin( patch->isZmin() )
{

    initElectroMagn3DQuantities( params, patch );
    
    // Charge currents currents and density for each species
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) { // end loop on ispec
        if ( emFields->Jx_s[ispec] != NULL ) {
            if ( emFields->Jx_s[ispec]->data_ != NULL )
                Jx_s[ispec]  = FieldFactory::create(dimPrim, 0, false, emFields->Jx_s[ispec]->name, params);
            else
                Jx_s[ispec]  = FieldFactory::create(emFields->Jx_s[ispec]->name, dimPrim, params);
        }
        if ( emFields->Jy_s[ispec] != NULL ) {
            if ( emFields->Jy_s[ispec]->data_ != NULL )
                Jy_s[ispec]  = FieldFactory::create(dimPrim, 1, false, emFields->Jy_s[ispec]->name, params);
            else
                Jy_s[ispec]  = FieldFactory::create(emFields->Jy_s[ispec]->name, dimPrim, params);
        }
        if ( emFields->Jz_s[ispec] != NULL ) {
            if ( emFields->Jz_s[ispec]->data_ != NULL )
                Jz_s[ispec]  = FieldFactory::create(dimPrim, 2, false, emFields->Jz_s[ispec]->name, params);
            else
                Jz_s[ispec]  = FieldFactory::create(emFields->Jz_s[ispec]->name, dimPrim, params);
        }
        if( emFields->rho_s[ispec] != NULL ) {
            if( emFields->rho_s[ispec]->data_ != NULL ) {
                rho_s[ispec] = new Field3D( dimPrim, emFields->rho_s[ispec]->name );
            } else {
                rho_s[ispec]  = new Field3D( emFields->rho_s[ispec]->name, dimPrim );
            }
        }
        
        if( params.Laser_Envelope_model ) {
            if( emFields->Env_Chi_s[ispec] != NULL ) {
                if( emFields->Env_Chi_s[ispec]->data_ != NULL ) {
                    Env_Chi_s[ispec] = new Field3D( dimPrim, emFields->Env_Chi_s[ispec]->name );
                } else {
                    Env_Chi_s[ispec]  = new Field3D( emFields->Env_Chi_s[ispec]->name, dimPrim );
                }
            }
        }
        
        
        
    } // loop on ispec
    
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize quantities used in ElectroMagn3D
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::initElectroMagn3DQuantities( Params &params, Patch *patch )
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
    for( size_t i=0 ; i<nDim_field ; i++ ) {
        // Standard scheme
        dimPrim[i] = n_space[i]+1;
        dimDual[i] = n_space[i]+2-(params.is_pxr);
        // + Ghost domain
        dimPrim[i] += 2*oversize[i];
        dimDual[i] += 2*oversize[i];
    }
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nx_d = n_space[0]+2+2*oversize[0]-(params.is_pxr);
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = n_space[1]+1+2*oversize[1];
    ny_d = n_space[1]+2+2*oversize[1]-(params.is_pxr);
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = n_space[2]+1+2*oversize[2];
    nz_d = n_space[2]+2+2*oversize[2]-(params.is_pxr);
    
    // Allocation of the EM fields
    Ex_  = FieldFactory::create( dimPrim, 0, false, "Ex", params );
    Ey_  = FieldFactory::create( dimPrim, 1, false, "Ey", params );
    Ez_  = FieldFactory::create( dimPrim, 2, false, "Ez", params );
    Bx_  = FieldFactory::create( dimPrim, 0, true,  "Bx", params );
    By_  = FieldFactory::create( dimPrim, 1, true,  "By", params );
    Bz_  = FieldFactory::create( dimPrim, 2, true,  "Bz", params );
    Bx_m = FieldFactory::create( dimPrim, 0, true,  "Bx_m", params );
    By_m = FieldFactory::create( dimPrim, 1, true,  "By_m", params );
    Bz_m = FieldFactory::create( dimPrim, 2, true,  "Bz_m", params );
    if( params.Laser_Envelope_model ) {
        Env_A_abs_ = new Field3D( dimPrim, "Env_A_abs" );
        Env_Chi_   = new Field3D( dimPrim, "Env_Chi" );
        Env_E_abs_ = new Field3D( dimPrim, "Env_E_abs" );
	Env_Ex_abs_= new Field3D( dimPrim, "Env_Ex_abs" );
    }
    
    // Total charge currents and densities
    Jx_   = FieldFactory::create( dimPrim, 0, false, "Jx", params );
    Jy_   = FieldFactory::create( dimPrim, 1, false, "Jy", params );
    Jz_   = FieldFactory::create( dimPrim, 2, false, "Jz", params );
    rho_  = new Field3D(dimPrim, "Rho" );
    
    //Edge coeffs are organized as follow and do not account for corner points
    //xmin/ymin - xmin/ymax - xmin/zmin - xmin/zmax - xmax/ymin - xmax/ymax - xmax/zmin - xmax/zmax
    //ymin/xmin - ymin/xmax - ymin/zmin - ymin/zmax - ymax/xmin - ymax/xmax - ymax/zmin - ymax/zmax
    //zmin/xmin - zmin/xmax - zmin/ymin - zmin/ymax - zmaz/xmin - zmaz/xmax - zmax/ymin - zmax/ymax
    beta_edge.resize( 24 );
    S_edge.resize( 24 );
    
        if(params.is_pxr == true) {
        rhoold_ = new Field3D( dimPrim, "RhoOld" );
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
     MESSAGE("index_bc_min / index_bc_max / nx_p / nx_d" << index_bc_min[0]
     << " " << index_bc_max[0] << " " << nx_p<< " " << nx_d);
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


void ElectroMagn3D::initPoisson( Patch *patch )
{
    Field3D *rho3D = static_cast<Field3D *>( rho_ );
    
    // Min and max indices for calculation of the scalar product (for primal & dual grid)
    //     scalar products are computed accounting only on real nodes
    //     ghost cells are used only for the (non-periodic) boundaries
    // dual indexes suppressed during "patchization"
    // ----------------------------------------------------------------------------------
    
    index_min_p_.resize( 3, 0 );
    index_max_p_.resize( 3, 0 );
    
    index_min_p_[0] = oversize[0];
    index_min_p_[1] = oversize[1];
    index_min_p_[2] = oversize[2];
    index_max_p_[0] = nx_p - 2 - oversize[0];
    index_max_p_[1] = ny_p - 2 - oversize[1];
    index_max_p_[2] = nz_p - 2 - oversize[2];
    if( patch->isXmin() ) {
        index_min_p_[0] = 0;
    }
    if( patch->isXmax() ) {
        index_max_p_[0] = nx_p-1;
    }
    
    phi_ = new Field3D( dimPrim );  // scalar potential
    r_   = new Field3D( dimPrim );  // residual vector
    p_   = new Field3D( dimPrim );  // direction vector
    Ap_  = new Field3D( dimPrim );  // A*p vector
    
    
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *phi_ )( i, j, k )   = 0.0;
                ( *r_ )( i, j, k )     = -( *rho3D )( i, j, k );
                ( *p_ )( i, j, k )     = ( *r_ )( i, j, k );
            }
        }//j
    }//i
    
} // initPoisson

double ElectroMagn3D::compute_r()
{
    double rnew_dot_rnew_local( 0. );
    for( unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++ ) {
        for( unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++ ) {
            for( unsigned int k=index_min_p_[2]; k<=index_max_p_[2]; k++ ) {
                rnew_dot_rnew_local += ( *r_ )( i, j, k )*( *r_ )( i, j, k );
            }
        }
    }
    return rnew_dot_rnew_local;
} // compute_r

void ElectroMagn3D::compute_Ap( Patch *patch )
{
    double one_ov_dx_sq       = 1.0/( dx*dx );
    double one_ov_dy_sq       = 1.0/( dy*dy );
    double one_ov_dz_sq       = 1.0/( dz*dz );
    double two_ov_dx2dy2dz2 = 2.0*( 1.0/( dx*dx )+1.0/( dy*dy )+1.0/( dz*dz ) );
    
    // vector product Ap = A*p
    for( unsigned int i=1; i<nx_p-1; i++ ) {
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( i, j, k ) = one_ov_dx_sq*( ( *p_ )( i-1, j, k )+( *p_ )( i+1, j, k ) )
                                      + one_ov_dy_sq*( ( *p_ )( i, j-1, k )+( *p_ )( i, j+1, k ) )
                                      + one_ov_dz_sq*( ( *p_ )( i, j, k-1 )+( *p_ )( i, j, k+1 ) )
                                      - two_ov_dx2dy2dz2*( *p_ )( i, j, k );
            }//k
        }//j
    }//i
    
    
    // Xmin BC
    if( patch->isXmin() ) {
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( 0, j, k )      = one_ov_dx_sq*( ( *p_ )( 1, j, k ) )
                                           +              one_ov_dy_sq*( ( *p_ )( 0, j-1, k )+( *p_ )( 0, j+1, k ) )
                                           +              one_ov_dz_sq*( ( *p_ )( 0, j, k-1 )+( *p_ )( 0, j, k+1 ) )
                                           -              two_ov_dx2dy2dz2*( *p_ )( 0, j, k );
            }
        }
        // at corners
        ( *Ap_ )( 0, 0, 0 )           = one_ov_dx_sq*( ( *p_ )( 1, 0, 0 ) ) // Xmin/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, 1, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, 0, 1 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( 0, 0, 0 );
        ( *Ap_ )( 0, ny_p-1, 0 )      = one_ov_dx_sq*( ( *p_ )( 1, ny_p-1, 0 ) ) // Xmin/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, ny_p-2, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, ny_p-1, 1 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( 0, ny_p-1, 0 );
        ( *Ap_ )( 0, 0, nz_p-1 )      = one_ov_dx_sq*( ( *p_ )( 1, 0, nz_p-1 ) ) // Xmin/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, 1, nz_p-1 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, 0, nz_p-2 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( 0, 0, nz_p-1 );
        ( *Ap_ )( 0, ny_p-1, nz_p-1 ) = one_ov_dx_sq*( ( *p_ )( 1, ny_p-1, nz_p-1 ) ) // Xmin/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, ny_p-2, nz_p-1 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, ny_p-1, nz_p-2 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( 0, ny_p-1, nz_p-1 );
    }
    
    // Xmax BC
    if( patch->isXmax() ) {
    
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( nx_p-1, j, k ) = one_ov_dx_sq*( ( *p_ )( nx_p-2, j, k ) )
                                           +              one_ov_dy_sq*( ( *p_ )( nx_p-1, j-1, k )+( *p_ )( nx_p-1, j+1, k ) )
                                           +              one_ov_dz_sq*( ( *p_ )( nx_p-1, j, k-1 )+( *p_ )( nx_p-1, j, k+1 ) )
                                           -              two_ov_dx2dy2dz2*( *p_ )( nx_p-1, j, k );
            }
        }
        // at corners
        ( *Ap_ )( nx_p-1, 0, 0 )      = one_ov_dx_sq*( ( *p_ )( nx_p-2, 0, 0 ) ) // Xmax/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, 1, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, 0, 1 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( nx_p-1, 0, 0 );
        ( *Ap_ )( nx_p-1, ny_p-1, 0 ) = one_ov_dx_sq*( ( *p_ )( nx_p-2, ny_p-1, 0 ) ) // Xmax/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, ny_p-2, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, ny_p-1, 1 ) )
                                        -                   two_ov_dx2dy2dz2*( *p_ )( nx_p-1, ny_p-1, 0 );
        ( *Ap_ )( nx_p-1, 0, nz_p-1 )      = one_ov_dx_sq*( ( *p_ )( nx_p-2, 0, 0 ) ) // Xmax/Ymin/Zmax
                                             +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, 1, nz_p-1 ) )
                                             +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, 0, nz_p-2 ) )
                                             -                   two_ov_dx2dy2dz2*( *p_ )( nx_p-1, 0, nz_p-1 );
        ( *Ap_ )( nx_p-1, ny_p-1, nz_p-1 ) = one_ov_dx_sq*( ( *p_ )( nx_p-2, ny_p-1, nz_p-1 ) ) // Xmax/Ymax/Zmax
                                             +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, ny_p-2, nz_p-1 ) )
                                             +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, ny_p-1, nz_p-2 ) )
                                             -                   two_ov_dx2dy2dz2*( *p_ )( nx_p-1, ny_p-1, nz_p-1 );
    }
    
} // compute_Ap

void ElectroMagn3D::compute_Ap_relativistic_Poisson( Patch *patch, double gamma_mean )
{

    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    double one_ov_dx_sq_ov_gamma_sq       = 1.0/( dx*dx )/( gamma_mean*gamma_mean );
    double one_ov_dy_sq                   = 1.0/( dy*dy );
    double one_ov_dz_sq                   = 1.0/( dz*dz );
    double two_ov_dxgam2dy2dz2            = 2.0*( 1.0/( dx*dx )/( gamma_mean*gamma_mean )+1.0/( dy*dy )+1.0/( dz*dz ) );
    
    // vector product Ap = A*p
    for( unsigned int i=1; i<nx_p-1; i++ ) {
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( i, j, k ) = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( i-1, j, k )+( *p_ )( i+1, j, k ) )
                                      + one_ov_dy_sq*( ( *p_ )( i, j-1, k )+( *p_ )( i, j+1, k ) )
                                      + one_ov_dz_sq*( ( *p_ )( i, j, k-1 )+( *p_ )( i, j, k+1 ) )
                                      - two_ov_dxgam2dy2dz2*( *p_ )( i, j, k );
            }//k
        }//j
    }//i
    
    
    // Xmin BC
    if( patch->isXmin() ) {
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( 0, j, k )      = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( 1, j, k ) )
                                           +              one_ov_dy_sq*( ( *p_ )( 0, j-1, k )+( *p_ )( 0, j+1, k ) )
                                           +              one_ov_dz_sq*( ( *p_ )( 0, j, k-1 )+( *p_ )( 0, j, k+1 ) )
                                           -              two_ov_dxgam2dy2dz2*( *p_ )( 0, j, k );
            }
        }
        // at corners
        ( *Ap_ )( 0, 0, 0 )           = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( 1, 0, 0 ) ) // Xmin/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, 1, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, 0, 1 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( 0, 0, 0 );
        ( *Ap_ )( 0, ny_p-1, 0 )      = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( 1, ny_p-1, 0 ) ) // Xmin/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, ny_p-2, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, ny_p-1, 1 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( 0, ny_p-1, 0 );
        ( *Ap_ )( 0, 0, nz_p-1 )      = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( 1, 0, nz_p-1 ) ) // Xmin/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, 1, nz_p-1 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, 0, nz_p-2 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( 0, 0, nz_p-1 );
        ( *Ap_ )( 0, ny_p-1, nz_p-1 ) = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( 1, ny_p-1, nz_p-1 ) ) // Xmin/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( 0, ny_p-2, nz_p-1 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( 0, ny_p-1, nz_p-2 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( 0, ny_p-1, nz_p-1 );
    }
    
    // Xmax BC
    if( patch->isXmax() ) {
    
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            for( unsigned int k=1; k<nz_p-1; k++ ) {
                ( *Ap_ )( nx_p-1, j, k ) = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( nx_p-2, j, k ) )
                                           +              one_ov_dy_sq*( ( *p_ )( nx_p-1, j-1, k )+( *p_ )( nx_p-1, j+1, k ) )
                                           +              one_ov_dz_sq*( ( *p_ )( nx_p-1, j, k-1 )+( *p_ )( nx_p-1, j, k+1 ) )
                                           -              two_ov_dxgam2dy2dz2*( *p_ )( nx_p-1, j, k );
            }
        }
        // at corners
        ( *Ap_ )( nx_p-1, 0, 0 )      = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( nx_p-2, 0, 0 ) ) // Xmax/Ymin/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, 1, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, 0, 1 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( nx_p-1, 0, 0 );
        ( *Ap_ )( nx_p-1, ny_p-1, 0 ) = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( nx_p-2, ny_p-1, 0 ) ) // Xmax/Ymax/Zmin
                                        +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, ny_p-2, 0 ) )
                                        +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, ny_p-1, 1 ) )
                                        -                   two_ov_dxgam2dy2dz2*( *p_ )( nx_p-1, ny_p-1, 0 );
        ( *Ap_ )( nx_p-1, 0, nz_p-1 )      = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( nx_p-2, 0, 0 ) ) // Xmax/Ymin/Zmax
                                             +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, 1, nz_p-1 ) )
                                             +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, 0, nz_p-2 ) )
                                             -                   two_ov_dxgam2dy2dz2*( *p_ )( nx_p-1, 0, nz_p-1 );
        ( *Ap_ )( nx_p-1, ny_p-1, nz_p-1 ) = one_ov_dx_sq_ov_gamma_sq*( ( *p_ )( nx_p-2, ny_p-1, nz_p-1 ) ) // Xmax/Ymax/Zmax
                                             +                   one_ov_dy_sq*( ( *p_ )( nx_p-1, ny_p-2, nz_p-1 ) )
                                             +                   one_ov_dz_sq*( ( *p_ )( nx_p-1, ny_p-1, nz_p-2 ) )
                                             -                   two_ov_dxgam2dy2dz2*( *p_ )( nx_p-1, ny_p-1, nz_p-1 );
    }
    
} // compute_Ap_relativistic_Poisson

double ElectroMagn3D::compute_pAp()
{
    double p_dot_Ap_local = 0.0;
    for( unsigned int i=index_min_p_[0]; i<=index_max_p_[0]; i++ ) {
        for( unsigned int j=index_min_p_[1]; j<=index_max_p_[1]; j++ ) {
            for( unsigned int k=index_min_p_[2]; k<=index_max_p_[2]; k++ ) {
                p_dot_Ap_local += ( *p_ )( i, j, k )*( *Ap_ )( i, j, k );
            }
        }
    }
    return p_dot_Ap_local;
} // compute_pAp

void ElectroMagn3D::update_pand_r( double r_dot_r, double p_dot_Ap )
{
    double alpha_k = r_dot_r/p_dot_Ap;
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *phi_ )( i, j, k ) += alpha_k * ( *p_ )( i, j, k );
                ( *r_ )( i, j, k )   -= alpha_k * ( *Ap_ )( i, j, k );
            }
        }
    }
    
} // update_pand_r

void ElectroMagn3D::update_p( double rnew_dot_rnew, double r_dot_r )
{
    double beta_k = rnew_dot_rnew/r_dot_r;
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *p_ )( i, j, k ) = ( *r_ )( i, j, k ) + beta_k * ( *p_ )( i, j, k );
            }
        }
    }
} // update_p

void ElectroMagn3D::initE( Patch *patch )
{
    Field3D *Ex3D  = static_cast<Field3D *>( Ex_ );
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_ );
    Field3D *rho3D = static_cast<Field3D *>( rho_ );
    
    // ------------------------------------------
    // Compute the electrostatic fields Ex and Ey
    // ------------------------------------------
    // Ex
    DEBUG( "Computing Ex from scalar potential" );
    for( unsigned int i=1; i<nx_d-1; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( i, j, k ) = ( ( *phi_ )( i-1, j, k )-( *phi_ )( i, j, k ) )/dx;
            }
        }
    }
    // Ey
    DEBUG( "Computing Ey from scalar potential" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=1; j<ny_d-1; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ey3D )( i, j, k ) = ( ( *phi_ )( i, j-1, k )-( *phi_ )( i, j, k ) )/dy;
            }
        }
    }
    DEBUG( "Computing Ez from scalar potential" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=1; k<nz_d-1; k++ ) {
                ( *Ez3D )( i, j, k ) = ( ( *phi_ )( i, j, k-1 )-( *phi_ )( i, j, k ) )/dz;
            }
        }
    }
    
    // Apply BC on Ex and Ey
    // ---------------------
    // Ex / Xmin
    if( patch->isXmin() ) {
        DEBUG( "Computing Xmin BC on Ex" );
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( 0, j, k ) = ( *Ex3D )( 1, j, k ) - dx*( *rho3D )( 0, j, k )
                                       + ( ( *Ey3D )( 0, j+1, k )-( *Ey3D )( 0, j, k ) )*dx/dy
                                       + ( ( *Ez3D )( 0, j, k+1 )-( *Ez3D )( 0, j, k ) )*dx/dz;
            }
        }
    }
    // Ex / Xmax
    if( patch->isXmax() ) {
        DEBUG( "Computing Xmax BC on Ex" );
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( nx_d-1, j, k ) = ( *Ex3D )( nx_d-2, j, k ) + dx*( *rho3D )( nx_p-1, j, k )
                                            - ( ( *Ey3D )( nx_p-1, j+1, k )-( *Ey3D )( nx_p-1, j, k ) )*dx/dy
                                            - ( ( *Ez3D )( nx_p-1, j, k+1 )-( *Ez3D )( nx_p-1, j, k ) )*dx/dz;
            }
        }
    }
    
    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;
    
} // initE

void ElectroMagn3D::initE_relativistic_Poisson( Patch *patch, double gamma_mean )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    Field3D *Ex3D  = static_cast<Field3D *>( Ex_rel_ );
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_rel_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_rel_ );
    Field3D *rho3D = static_cast<Field3D *>( rho_ );
    
    // ------------------------------------------
    // Compute the fields Ex, Ey and Ez
    // ------------------------------------------
    
    
    
    // Ex
    MESSAGE( 1, "Computing Ex from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=1; i<nx_d-1; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( i, j, k ) = ( ( *phi_ )( i-1, j, k )-( *phi_ )( i, j, k ) )/dx/gamma_mean/gamma_mean;
            }
        }
    }
    MESSAGE( 1, "Ex: done" );
    // Ey
    MESSAGE( 1, "Computing Ey from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=1; j<ny_d-1; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ey3D )( i, j, k ) = ( ( *phi_ )( i, j-1, k )-( *phi_ )( i, j, k ) )/dy;
            }
        }
    }
    MESSAGE( 1, "Ey: done" );
    // Ez
    MESSAGE( 1, "Computing Ez from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=1; k<nz_d-1; k++ ) {
                ( *Ez3D )( i, j, k ) = ( ( *phi_ )( i, j, k-1 )-( *phi_ )( i, j, k ) )/dz;
            }
        }
    }
    MESSAGE( 1, "Ez: done" );
    
    // Apply BC on Ex, Ey and Ez
    // ---------------------
    // Ex / Xmin
    if( patch->isXmin() ) {
        DEBUG( "Computing Xmin BC on Ex, relativistic Poisson problem" );
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( 0, j, k ) = ( *Ex3D )( 1, j, k ) - dx*( *rho3D )( 0, j, k )
                                       + ( ( *Ey3D )( 0, j+1, k )-( *Ey3D )( 0, j, k ) )*dx/dy
                                       + ( ( *Ez3D )( 0, j, k+1 )-( *Ez3D )( 0, j, k ) )*dx/dz;
            }
        }
    }
    // Ex / Xmax
    if( patch->isXmax() ) {
        DEBUG( "Computing Xmax BC on Ex, relativistic Poisson problem" );
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( nx_d-1, j, k ) = ( *Ex3D )( nx_d-2, j, k ) + dx*( *rho3D )( nx_p-1, j, k )
                                            - ( ( *Ey3D )( nx_p-1, j+1, k )-( *Ey3D )( nx_p-1, j, k ) )*dx/dy
                                            - ( ( *Ez3D )( nx_p-1, j, k+1 )-( *Ez3D )( nx_p-1, j, k ) )*dx/dz;
            }
        }
    }
    
    // // Ey / Ymin
    // if (patch->isYmin()) {
    //     DEBUG("Computing Ymin BC on Ey, relativistic Poisson problem");
    //     for (unsigned int i=0; i<nx_p; i++) {
    //         for (unsigned int k=0; k<nz_p; k++) {
    //             (*Ey3D)(i,0,k) = (*Ey3D)(i,1,k) - dy*(*rho3D)(i,0,k)
    //                 + ((*Ex3D)(i+1,0,k)-(*Ex3D)(i,0,k))*dy/dx
    //                 + ((*Ez3D)(i,0,k+1)-(*Ez3D)(i,0,k))*dy/dz;
    //         }
    //     }
    // }
    //
    // // Ey / Ymax
    // if (patch->isYmax()) {
    //     DEBUG("Computing Ymax BC on Ey, relativistic Poisson problem");
    //     for (unsigned int i=0; i<nx_p; i++) {
    //         for (unsigned int k=0; k<nz_p; k++) {
    //             (*Ey3D)(i,ny_d-1,k) = (*Ey3D)(i,ny_d-2,k) + dy*(*rho3D)(i,ny_p-1,k)
    //                 - ((*Ex3D)(i+1,ny_p-1,k)-(*Ex3D)(i,ny_p-1,k))*dy/dx
    //                 - ((*Ez3D)(i,ny_p-1,k+1)-(*Ez3D)(i,ny_p-1,k))*dy/dz;
    //         }
    //     }
    // }
    //
    // // Ez / Zmin
    // if (patch->isZmin()) {
    //     DEBUG("Computing Zmin BC on Ez, relativistic Poisson problem");
    //     for (unsigned int i=0; i<nx_p; i++) {
    //         for (unsigned int j=0; j<ny_p; j++) {
    //             (*Ez3D)(i,j,0) = (*Ez3D)(i,j,1) - dz*(*rho3D)(i,j,0)
    //                 + ((*Ey3D)(i,j+1,0)-(*Ey3D)(i,j,0))*dz/dy
    //                 + ((*Ex3D)(i+1,j,0)-(*Ex3D)(i,j,0))*dz/dx;
    //         }
    //     }
    // }
    //
    // // Ez / Zmax
    // if (patch->isZmax()) {
    //     DEBUG("Computing Zmax BC on Ez, relativistic Poisson problem");
    //     for (unsigned int i=0; i<nx_p; i++) {
    //         for (unsigned int j=0; j<ny_p; j++) {
    //             (*Ez3D)(i,j,nz_d-1) = (*Ez3D)(i,j,nz_d-2) + dz*(*rho3D)(i,j,nz_p-1)
    //                 - ((*Ey3D)(i,j+1,nz_p-1)-(*Ey3D)(i,j,nz_p-1))*dz/dy
    //                 - ((*Ex3D)(i+1,j,nz_p-1)-(*Ex3D)(i,j,nz_p-1))*dz/dx;
    //         }
    //     }
    // }
    
    delete phi_;
    delete r_;
    delete p_;
    delete Ap_;
    
} // initE_relativistic_Poisson

void ElectroMagn3D::initB_relativistic_Poisson( Patch *patch, double gamma_mean )
{
    // gamma_mean is the average Lorentz factor of the species whose fields will be computed
    // See for example https://doi.org/10.1016/j.nima.2016.02.043 for more details
    
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_rel_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_rel_ );
    
    // Bx is zero everywhere
    Field3D *Bx3D  = static_cast<Field3D *>( Bx_rel_ );
    Field3D *By3D  = static_cast<Field3D *>( By_rel_ );
    Field3D *Bz3D  = static_cast<Field3D *>( Bz_rel_ );
    
    double beta_mean = sqrt( 1.-1./gamma_mean/gamma_mean );
    
    // ------------------------------------------
    // Compute the fields Ex and Ey
    // ------------------------------------------
    // Bx is identically zero
    // (hypothesis of negligible J transverse with respect to Jx)
    MESSAGE( 1, "Computing Bx relativistic Poisson problem" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *Bx3D )( i, j, k ) = 0.;
            }
        }
    }
    MESSAGE( 1, "Bx: done" );
    
    // By
    MESSAGE( 1, "Computing By from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *By3D )( i, j, k ) = -beta_mean*( *Ez3D )( i, j, k );
            }
        }
    }
    MESSAGE( 1, "By: done" );
    // Bz
    MESSAGE( 1, "Computing Bz from scalar potential, relativistic Poisson problem" );
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Bz3D )( i, j, k ) = beta_mean*( *Ey3D )( i, j, k );
            }
        }
    }
    MESSAGE( 1, "Bz: done" );
    
    
    
} // initB_relativistic_Poisson

void ElectroMagn3D::initRelativisticPoissonFields( Patch *patch )
{
    // ------ Init temporary fields for relativistic field initialization
    
    // E fields centered as in FDTD, to be added to the already present electric fields
    Ex_rel_  = new Field3D( dimPrim, 0, false, "Ex_rel" );
    Ey_rel_  = new Field3D( dimPrim, 1, false, "Ey_rel" );
    Ez_rel_  = new Field3D( dimPrim, 2, false, "Ez_rel" );
    
    
    // B fields centered as the E fields in FDTD (Bx null)
    Bx_rel_  = new Field3D( dimPrim, 0, true,  "Bx_rel" ); // will be identically zero (hypothesis of negligible transverse current with respect to longitudinal current)
    By_rel_  = new Field3D( dimPrim, 2, false,  "By_rel" ); // is equal to -beta*Ez, thus it inherits the same centering of Ez
    Bz_rel_  = new Field3D( dimPrim, 1, false,  "Bz_rel" ); // is equal to  beta*Ey, thus it inherits the same centering of Ey
    
    
    // ----- B fields centered as in FDTD, to be added to the already present magnetic fields
    
    // B field advanced by dt/2
    Bx_rel_t_plus_halfdt_  = new Field3D( dimPrim, 0, true,  "Bx_rel_t_plus_halfdt" );
    By_rel_t_plus_halfdt_  = new Field3D( dimPrim, 1, true,  "By_rel_t_plus_halfdt" );
    Bz_rel_t_plus_halfdt_  = new Field3D( dimPrim, 2, true,  "Bz_rel_t_plus_halfdt" );
    // B field "advanced" by -dt/2
    Bx_rel_t_minus_halfdt_  = new Field3D( dimPrim, 0, true,  "Bx_rel_t_plus_halfdt" );
    By_rel_t_minus_halfdt_  = new Field3D( dimPrim, 1, true,  "By_rel_t_plus_halfdt" );
    Bz_rel_t_minus_halfdt_  = new Field3D( dimPrim, 2, true,  "Bz_rel_t_plus_halfdt" );
    
    
} // initRelativisticPoissonFields

void ElectroMagn3D::sum_rel_fields_to_em_fields( Patch *patch )
{
    Field3D *Ex3Drel  = static_cast<Field3D *>( Ex_rel_ );
    Field3D *Ey3Drel  = static_cast<Field3D *>( Ey_rel_ );
    Field3D *Ez3Drel  = static_cast<Field3D *>( Ez_rel_ );
    
    // B_t_plus_halfdt
    Field3D *Bx_rel_t_plus_halfdt = static_cast<Field3D *>( Bx_rel_t_plus_halfdt_ );
    Field3D *By_rel_t_plus_halfdt = static_cast<Field3D *>( By_rel_t_plus_halfdt_ );
    Field3D *Bz_rel_t_plus_halfdt = static_cast<Field3D *>( Bz_rel_t_plus_halfdt_ );
    
    // B_t_minus_halfdt
    Field3D *Bx_rel_t_minus_halfdt = static_cast<Field3D *>( Bx_rel_t_minus_halfdt_ );
    Field3D *By_rel_t_minus_halfdt = static_cast<Field3D *>( By_rel_t_minus_halfdt_ );
    Field3D *Bz_rel_t_minus_halfdt = static_cast<Field3D *>( Bz_rel_t_minus_halfdt_ );
    
    // E and B fields already existing on the grid
    Field3D *Ex3D  = static_cast<Field3D *>( Ex_ );
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_ );
    Field3D *Bx3D  = static_cast<Field3D *>( Bx_ );
    Field3D *By3D  = static_cast<Field3D *>( By_ );
    Field3D *Bz3D  = static_cast<Field3D *>( Bz_ );
    Field3D *Bx3D0  = static_cast<Field3D *>( Bx_m );
    Field3D *By3D0  = static_cast<Field3D *>( By_m );
    Field3D *Bz3D0  = static_cast<Field3D *>( Bz_m );
    
    // Ex (d,p,p)
    for( unsigned int i=0; i<nx_d; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( i, j, k ) = ( *Ex3D )( i, j, k ) + ( *Ex3Drel )( i, j, k );
            }
        }
    }
    
    // Ey (p,d,p)
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ey3D )( i, j, k ) = ( *Ey3D )( i, j, k ) + ( *Ey3Drel )( i, j, k );
            }
        }
    }
    
    // Ez (p,p,d)
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *Ez3D )( i, j, k ) = ( *Ez3D )( i, j, k ) + ( *Ez3Drel )( i, j, k );
            }
        }
    }
    
    
    
    // Since Brel is centered in time as E, it is inconsistent with FDTD,
    // where E and B are staggered in time.
    // Possible solution:
    // Use FDTD scheme to integrate Maxwell-Faraday equation forward in time by dt/2 to obtain B
    // Use FDTD scheme to integrate Maxwell-Faraday equation backwards in time by dt/2 to obtain Bm
    // Add the forward-evolved and backward-evolved fields to the grid fields
    
    double half_dt_ov_dx = 0.5 * timestep / dx;
    double half_dt_ov_dy = 0.5 * timestep / dy;
    double half_dt_ov_dz = 0.5 * timestep / dz;
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                // forward advance by dt/2
                ( *Bx_rel_t_plus_halfdt )( i, j, k ) += -half_dt_ov_dy * ( ( *Ez3Drel )( i, j, k ) - ( *Ez3Drel )( i, j-1, k ) ) + half_dt_ov_dz * ( ( *Ey3Drel )( i, j, k ) - ( *Ey3Drel )( i, j, k-1 ) );
                // backward advance by dt/2
                ( *Bx_rel_t_minus_halfdt )( i, j, k ) -= -half_dt_ov_dy * ( ( *Ez3Drel )( i, j, k ) - ( *Ez3Drel )( i, j-1, k ) ) + half_dt_ov_dz * ( ( *Ey3Drel )( i, j, k ) - ( *Ey3Drel )( i, j, k-1 ) );
                // sum to the fields on grid
                ( *Bx3D )( i, j, k ) += ( *Bx_rel_t_plus_halfdt )( i, j, k );
                ( *Bx3D0 )( i, j, k ) += ( *Bx_rel_t_minus_halfdt )( i, j, k );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                // forward advance by dt/2
                ( *By_rel_t_plus_halfdt )( i, j, k ) += -half_dt_ov_dz * ( ( *Ex3Drel )( i, j, k ) - ( *Ex3Drel )( i, j, k-1 ) ) + half_dt_ov_dx * ( ( *Ez3Drel )( i, j, k ) - ( *Ez3Drel )( i-1, j, k ) );
                // backward advance by dt/2
                ( *By_rel_t_minus_halfdt )( i, j, k ) -= -half_dt_ov_dz * ( ( *Ex3Drel )( i, j, k ) - ( *Ex3Drel )( i, j, k-1 ) ) + half_dt_ov_dx * ( ( *Ez3Drel )( i, j, k ) - ( *Ez3Drel )( i-1, j, k ) );
                // sum to the fields on grid
                ( *By3D )( i, j, k ) += ( *By_rel_t_plus_halfdt )( i, j, k );
                ( *By3D0 )( i, j, k ) += ( *By_rel_t_minus_halfdt )( i, j, k );
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                // forward advance by dt/2
                ( *Bz_rel_t_plus_halfdt )( i, j, k ) += -half_dt_ov_dx * ( ( *Ey3Drel )( i, j, k ) - ( *Ey3Drel )( i-1, j, k ) ) + half_dt_ov_dy * ( ( *Ex3Drel )( i, j, k ) - ( *Ex3Drel )( i, j-1, k ) );
                // backward advance by dt/2
                ( *Bz_rel_t_minus_halfdt )( i, j, k ) -= -half_dt_ov_dx * ( ( *Ey3Drel )( i, j, k ) - ( *Ey3Drel )( i-1, j, k ) ) + half_dt_ov_dy * ( ( *Ex3Drel )( i, j, k ) - ( *Ex3Drel )( i, j-1, k ) );
                // sum to the fields on grid
                ( *Bz3D )( i, j, k ) += ( *Bz_rel_t_plus_halfdt )( i, j, k );
                ( *Bz3D0 )( i, j, k ) += ( *Bz_rel_t_minus_halfdt )( i, j, k );
            }
        }
    }
    
    
    // delete temporary fields used for relativistic initialization
    delete Ex_rel_;
    delete Ey_rel_;
    delete Ez_rel_;
    delete Bx_rel_;
    delete By_rel_;
    delete Bz_rel_;
    
    delete Bx_rel_t_plus_halfdt;
    delete By_rel_t_plus_halfdt;
    delete Bz_rel_t_plus_halfdt;
    delete Bx_rel_t_minus_halfdt;
    delete By_rel_t_minus_halfdt;
    delete Bz_rel_t_minus_halfdt;
    
    
} // sum_rel_fields_to_em_fields


void ElectroMagn3D::centeringE( std::vector<double> E_Add )
{
    Field3D *Ex3D  = static_cast<Field3D *>( Ex_ );
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_ );
    
    // Centering electrostatic fields
    for( unsigned int i=0; i<nx_d; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( i, j, k ) += E_Add[0];
            }
        }
    }
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ey3D )( i, j, k ) += E_Add[1];
            }
        }
    }
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *Ez3D )( i, j, k ) += E_Add[2];
            }
        }
    }
    
} // centeringE

void ElectroMagn3D::centeringErel( std::vector<double> E_Add )
{
    Field3D *Ex3D  = static_cast<Field3D *>( Ex_rel_ );
    Field3D *Ey3D  = static_cast<Field3D *>( Ey_rel_ );
    Field3D *Ez3D  = static_cast<Field3D *>( Ez_rel_ );
    
    // Centering electrostatic fields
    for( unsigned int i=0; i<nx_d; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ex3D )( i, j, k ) += E_Add[0];
            }
        }
    }
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Ey3D )( i, j, k ) += E_Add[1];
            }
        }
    }
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *Ez3D )( i, j, k ) += E_Add[2];
            }
        }
    }
    
} // centeringErel

// ---------------------------------------------------------------------------------------------------------------------
// End of Solve Poisson methods
// ---------------------------------------------------------------------------------------------------------------------


// ---------------------------------------------------------------------------------------------------------------------
// Save the former Magnetic-Fields (used to center them)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::saveMagneticFields( bool is_spectral )
{
    // Static cast of the fields
    if( !is_spectral ) {
        Field3D *Bx3D   = static_cast<Field3D *>( Bx_ );
        Field3D *By3D   = static_cast<Field3D *>( By_ );
        Field3D *Bz3D   = static_cast<Field3D *>( Bz_ );
        Field3D *Bx3D_m = static_cast<Field3D *>( Bx_m );
        Field3D *By3D_m = static_cast<Field3D *>( By_m );
        Field3D *Bz3D_m = static_cast<Field3D *>( Bz_m );
        
        // Magnetic field Bx^(p,d,d)
        memcpy( &( ( *Bx3D_m )( 0, 0, 0 ) ), &( ( *Bx3D )( 0, 0, 0 ) ), nx_p*ny_d*nz_d*sizeof( double ) );
        
        // Magnetic field By^(d,p,d)
        memcpy( &( ( *By3D_m )( 0, 0, 0 ) ), &( ( *By3D )( 0, 0, 0 ) ), nx_d*ny_p*nz_d*sizeof( double ) );
        
        // Magnetic field Bz^(d,d,p)
        memcpy( &( ( *Bz3D_m )( 0, 0, 0 ) ), &( ( *Bz3D )( 0, 0, 0 ) ), nx_d*ny_d*nz_p*sizeof( double ) );
    } else {
        Bx_m->deallocateDataAndSetTo( Bx_ );
        By_m->deallocateDataAndSetTo( By_ );
        Bz_m->deallocateDataAndSetTo( Bz_ );
    }
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
Field *ElectroMagn3D::createField( string fieldname, Params& params )
{
    if     (fieldname.substr(0,2)=="Ex" ) return FieldFactory::create(dimPrim, 0, false, fieldname, params);
    else if(fieldname.substr(0,2)=="Ey" ) return FieldFactory::create(dimPrim, 1, false, fieldname, params);
    else if(fieldname.substr(0,2)=="Ez" ) return FieldFactory::create(dimPrim, 2, false, fieldname, params);
    else if(fieldname.substr(0,2)=="Bx" ) return FieldFactory::create(dimPrim, 0, true,  fieldname, params);
    else if(fieldname.substr(0,2)=="By" ) return FieldFactory::create(dimPrim, 1, true,  fieldname, params);
    else if(fieldname.substr(0,2)=="Bz" ) return FieldFactory::create(dimPrim, 2, true,  fieldname, params);
    else if(fieldname.substr(0,2)=="Jx" ) return FieldFactory::create(dimPrim, 0, false, fieldname, params);
    else if(fieldname.substr(0,2)=="Jy" ) return FieldFactory::create(dimPrim, 1, false, fieldname, params);
    else if(fieldname.substr(0,2)=="Jz" ) return FieldFactory::create(dimPrim, 2, false, fieldname, params);
    else if(fieldname.substr(0,3)=="Rho") return new Field3D(dimPrim, fieldname );
    else if(fieldname.substr(0,9)=="Env_A_abs" ) return new Field3D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,7)=="Env_Chi" ) return new Field3D(dimPrim, 0, false, fieldname);
    else if(fieldname.substr(0,9)=="Env_E_abs" ) return new Field3D(dimPrim, 0, false, fieldname);

    ERROR("Cannot create field "<<fieldname);
    return NULL;
}


// ---------------------------------------------------------------------------------------------------------------------
// Center the Magnetic Fields (used to push the particle)
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::centerMagneticFields()
{
    // Static cast of the fields
    Field3D *Bx3D   = static_cast<Field3D *>( Bx_ );
    Field3D *By3D   = static_cast<Field3D *>( By_ );
    Field3D *Bz3D   = static_cast<Field3D *>( Bz_ );
    Field3D *Bx3D_m = static_cast<Field3D *>( Bx_m );
    Field3D *By3D_m = static_cast<Field3D *>( By_m );
    Field3D *Bz3D_m = static_cast<Field3D *>( Bz_m );
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            for( unsigned int k=0 ; k<nz_d ; k++ ) {
                ( *Bx3D_m )( i, j, k ) = ( ( *Bx3D )( i, j, k ) + ( *Bx3D_m )( i, j, k ) )*0.5;
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_d ; k++ ) {
                ( *By3D_m )( i, j, k ) = ( ( *By3D )( i, j, k ) + ( *By3D_m )( i, j, k ) )*0.5;
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D_m )( i, j, k ) = ( ( *Bz3D )( i, j, k ) + ( *Bz3D_m )( i, j, k ) )*0.5;
            }
        } // end for j
    } // end for i
    
    
}//END centerMagneticFields


// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass binomial filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::binomialCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes)
{
    // Static-cast of the currents
    Field3D *Jx3D = static_cast<Field3D *>( Jx_ );
    Field3D *Jy3D = static_cast<Field3D *>( Jy_ );
    Field3D *Jz3D = static_cast<Field3D *>( Jz_ );
    
    // applying a single pass of the binomial filter along X
    if (ipass < passes[0]){
        // on Jx^(d,p) -- external points are treated by exchange. Boundary points not concerned by exchange are treated with a lower order filter.
        for( unsigned int i=0; i<nx_d-1; i++ ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i+1, j, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=nx_d-2; i>0; i-- ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i-1, j, k ) )*0.5;
                }
            }
        }
        // Jy
        for( unsigned int i=0; i<nx_p-1; i++ ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i+1, j, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=nx_p-2; i>0; i-- ) {
            for( unsigned int j=0; j<ny_d; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i-1, j, k ) )*0.5;
                }
            }
        }
        // Jz
        for( unsigned int i=0; i<nx_p-1; i++ ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=0; k<nz_d; k++ ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i+1, j, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=nx_p-2; i>0; i-- ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=0; k<nz_d; k++ ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i-1, j, k ) )*0.5;
                }
            }
        }
    }

    // applying a single pass of the binomial filter along Y
    if (ipass < passes[1]){
        //Jx
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=0; j<ny_p-1; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i, j+1, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=ny_p-2; j>0; j-- ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i, j-1, k ) )*0.5;
                }
            }
        }
        //Jy
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=0; j<ny_d-1; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i, j+1, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=ny_d-2; j>0; j-- ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i, j-1, k ) )*0.5;
                }
            }
        }
        //Jz
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=0; j<ny_p-1; j++ ) {
                for( unsigned int k=0; k<nz_d; k++ ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i, j+1, k ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=ny_p-2; j>0; j-- ) {
                for( unsigned int k=0; k<nz_d; k++ ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i, j-1, k ) )*0.5;
                }
            }
        }
    }
    // applying a single pass of the binomial filter along Z
    if (ipass < passes[2]){
        //Jx
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=0; k<nz_p-1; k++ ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i, j, k+1 ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=nz_p-2; k>0; k-- ) {
                    ( *Jx3D )( i, j, k ) = ( ( *Jx3D )( i, j, k ) + ( *Jx3D )( i, j, k-1 ) )*0.5;
                }
            }
        }
        //Jy
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_d-1; j++ ) {
                for( unsigned int k=0; k<nz_p-1; k++ ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i, j, k+1 ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_d-1; j++ ) {
                for( unsigned int k=nz_p-2; k>0; k-- ) {
                    ( *Jy3D )( i, j, k ) = ( ( *Jy3D )( i, j, k ) + ( *Jy3D )( i, j, k-1 ) )*0.5;
                }
            }
        }
        //Jz
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=0; k<nz_d-1; k++ ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i, j, k+1 ) )*0.5;
                }
            }
        }
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=nz_d-2; k>0; k-- ) {
                    ( *Jz3D )( i, j, k ) = ( ( *Jz3D )( i, j, k ) + ( *Jz3D )( i, j, k-1 ) )*0.5;
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply a single pass custom FIR based filter on currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::customFIRCurrentFilter(unsigned int ipass, std::vector<unsigned int> passes, std::vector<double> filtering_coeff)
{
    // Static-cast of the currents
    Field3D *Jx3D = static_cast<Field3D *>( Jx_ );
    Field3D *Jy3D = static_cast<Field3D *>( Jy_ );
    Field3D *Jz3D = static_cast<Field3D *>( Jz_ );

    unsigned int m=1 ;

    // Guard-Cell Current
    unsigned int gcfilt=0 ;

    // Applying a single pass of the custom FIR based filter along X
    if (ipass < passes[0]){
        Field3D *tmp   = new Field3D( dimPrim, 0, false );
        tmp->copyFrom( Jx3D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_d-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=1; k<nz_p-1; k++ ) {
                    ( *Jx3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jx3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m, j, k ) ;
                    }
                    ( *Jx3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        tmp   = new Field3D( dimPrim, 1, false );
        tmp->copyFrom( Jy3D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            for( unsigned int j=1; j<ny_d-1; j++ ) {
                for( unsigned int k=1; k<nz_p-1; k++ ) {
                    ( *Jy3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jy3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m, j, k ) ;
                    }
                    ( *Jy3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        tmp   = new Field3D( dimPrim, 2, false );
        tmp->copyFrom( Jz3D );
        for( unsigned int i=((filtering_coeff.size()-1)/(m*2)+gcfilt); i<nx_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=1; k<nz_d-1; k++ ) {
                    ( *Jz3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jz3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m, j, k ) ;
                    }
                    ( *Jz3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
    }

    // Applying a single pass of the custom FIR based filter along Y
    if (ipass < passes[1]){
        // On Jx^(d,p,p) -- External points are treated by exchange
        Field3D *tmp   = new Field3D( dimPrim, 0, false );
        tmp->copyFrom( Jx3D );
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=((filtering_coeff.size()-1)/(m*2)+gcfilt); j<ny_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); j++ ) {
                for( unsigned int k=1; k<nz_p-1; k++ ) {
                    ( *Jx3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jx3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m, k ) ;
                    }
                    ( *Jx3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        // On Jy^(p,d,p) -- External points are treated by exchange
        tmp   = new Field3D( dimPrim, 1, false );
        tmp->copyFrom( Jy3D );
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=((filtering_coeff.size()-1)/(m*2)+gcfilt); j<ny_d-((filtering_coeff.size()-1)/(m*2)+gcfilt); j++ ) {
                for( unsigned int k=1; k<nz_p-1; k++ ) {
                    ( *Jy3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jy3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m, k ) ;
                    }
                    ( *Jy3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        // On Jz^(p,p,d) -- External points are treated by exchange
        tmp   = new Field3D( dimPrim, 2, false );
        tmp->copyFrom( Jz3D );
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=((filtering_coeff.size()-1)/(m*2)+gcfilt); j<ny_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); j++ ) {
                for( unsigned int k=1; k<nz_d-1; k++ ) {
                    ( *Jz3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jz3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m , k) ;
                    }
                    ( *Jz3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
    }

    // Applying a single pass of the custom FIR based filter along Z
    if (ipass < passes[2]){
        // On Jx^(d,p,p) -- External points are treated by exchange
        Field3D *tmp   = new Field3D( dimPrim, 0, false );
        tmp->copyFrom( Jx3D );
        for( unsigned int i=1; i<nx_d-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=((filtering_coeff.size()-1)/(m*2)+gcfilt); k<ny_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); k++ ) {
                    ( *Jx3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jx3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j, k - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
                    }
                    ( *Jx3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        // On Jy^(p,d,p) -- External points are treated by exchange
        tmp   = new Field3D( dimPrim, 1, false );
        tmp->copyFrom( Jy3D );
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_d-1; j++ ) {
                for( unsigned int k=((filtering_coeff.size()-1)/(m*2)+gcfilt); k<ny_p-((filtering_coeff.size()-1)/(m*2)+gcfilt); k++ ) {
                    ( *Jy3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jy3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j, k - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
                    }
                    ( *Jy3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
        // On Jz^(p,p,d) -- External points are treated by exchange
        tmp   = new Field3D( dimPrim, 2, false );
        tmp->copyFrom( Jz3D );
        for( unsigned int i=1; i<nx_p-1; i++ ) {
            for( unsigned int j=1; j<ny_p-1; j++ ) {
                for( unsigned int k=((filtering_coeff.size()-1)/(m*2)+gcfilt); k<ny_d-((filtering_coeff.size()-1)/(m*2)+gcfilt); k++ ) {
                    ( *Jz3D )( i, j, k ) = 0. ;
                    for ( unsigned int kernel_idx = 0; kernel_idx < filtering_coeff.size(); kernel_idx+=m) {
                        ( *Jz3D )( i, j, k ) += filtering_coeff[kernel_idx]*( *tmp )( i, j, k - (filtering_coeff.size()-1)/(m*2) + kernel_idx/m ) ;
                    }
                    ( *Jz3D )( i, j, k ) *= m ;
                }
            }
        }
        delete tmp;
    }

}//END customFIRCurrentFilter

void ElectroMagn3D::center_fields_from_relativistic_Poisson( Patch *patch )
{


    // B field centered in time as E field, at time t
    Field3D *Bx3Drel  = static_cast<Field3D *>( Bx_rel_ );
    Field3D *By3Drel  = static_cast<Field3D *>( By_rel_ );
    Field3D *Bz3Drel  = static_cast<Field3D *>( Bz_rel_ );
    
    // B field centered in time at time t+dt/2
    Field3D *Bx3D  = static_cast<Field3D *>( Bx_rel_t_plus_halfdt_ );
    Field3D *By3D  = static_cast<Field3D *>( By_rel_t_plus_halfdt_ );
    Field3D *Bz3D  = static_cast<Field3D *>( Bz_rel_t_plus_halfdt_ );
    // B field centered in time at time t-dt/2
    Field3D *Bx3D0  = static_cast<Field3D *>( Bx_rel_t_minus_halfdt_ );
    Field3D *By3D0  = static_cast<Field3D *>( By_rel_t_minus_halfdt_ );
    Field3D *Bz3D0  = static_cast<Field3D *>( Bz_rel_t_minus_halfdt_ );
    
    
    // The B_rel fields, centered as B, will be advanced by dt/2 and -dt/2
    // for proper centering in FDTD, but first they have to be centered in space
    // The advance by dt and -dt and the sum to the existing grid fields is performed in
    // ElectroMagn3D::sum_rel_fields_to_em_fields
    
    // Bx (p,d,d)   Bx_rel is identically zero and centered as Bx, no special interpolation of indices
    for( unsigned int i=0; i<nx_p; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *Bx3D )( i, j, k )= ( *Bx3Drel )( i, j, k );
                ( *Bx3D0 )( i, j, k )= ( *Bx3Drel )( i, j, k );
            }
        }
    }
    
    // ---------- center the B fields
    // By (d,p,d) - remember that Byrel is centered as Ezrel (p,p,d)
    for( unsigned int i=1; i<nx_d-1; i++ ) {
        for( unsigned int j=0; j<ny_p; j++ ) {
            for( unsigned int k=0; k<nz_d; k++ ) {
                ( *By3D )( i, j, k )= 0.5 * ( ( *By3Drel )( i, j, k ) + ( *By3Drel )( i-1, j, k ) );
                ( *By3D0 )( i, j, k )= 0.5 * ( ( *By3Drel )( i, j, k ) + ( *By3Drel )( i-1, j, k ) );
            }
        }
    }
    
    // Bz (d,d,p) - remember that Bzrel is centered as Eyrel (p,d,p)
    for( unsigned int i=1; i<nx_d-1; i++ ) {
        for( unsigned int j=0; j<ny_d; j++ ) {
            for( unsigned int k=0; k<nz_p; k++ ) {
                ( *Bz3D )( i, j, k )= 0.5 * ( ( *Bz3Drel )( i, j, k ) + ( *Bz3Drel )( i-1, j, k ) );
                ( *Bz3D0 )( i, j, k )= 0.5 * ( ( *Bz3Drel )( i, j, k ) + ( *Bz3Drel )( i-1, j, k ) );
            }
        }
    }
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total density and currents from species density and currents
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computeTotalRhoJ()
{
    // static cast of the total currents and densities
    Field3D *Jx3D    = static_cast<Field3D *>( Jx_ );
    Field3D *Jy3D    = static_cast<Field3D *>( Jy_ );
    Field3D *Jz3D    = static_cast<Field3D *>( Jz_ );
    Field3D *rho3D   = static_cast<Field3D *>( rho_ );
    
    
    // -----------------------------------
    // Species currents and charge density
    // -----------------------------------
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Jx_s[ispec] ) {
            Field3D *Jx3D_s  = static_cast<Field3D *>( Jx_s[ispec] );
            for( unsigned int i=0 ; i<Jx3D->dims_[0] ; i++ )
                for( unsigned int j=0 ; j<Jx3D->dims_[1] ; j++ )
                    for( unsigned int k=0 ; k<Jx3D->dims_[2] ; k++ ) {
                        ( *Jx3D )( i, j, k ) += ( *Jx3D_s )( i, j, k );
                    }
        }
        if( Jy_s[ispec] ) {
            Field3D *Jy3D_s  = static_cast<Field3D *>( Jy_s[ispec] );
            for( unsigned int i=0 ; i<Jy3D->dims_[0] ; i++ )
                for( unsigned int j=0 ; j<Jy3D->dims_[1] ; j++ )
                    for( unsigned int k=0 ; k<Jy3D->dims_[2] ; k++ ) {
                        ( *Jy3D )( i, j, k ) += ( *Jy3D_s )( i, j, k );
                    }
        }
        if( Jz_s[ispec] ) {
            Field3D *Jz3D_s  = static_cast<Field3D *>( Jz_s[ispec] );
            for( unsigned int i=0 ; i<Jz3D->dims_[0] ; i++ )
                for( unsigned int j=0 ; j<Jz3D->dims_[1] ; j++ )
                    for( unsigned int k=0 ; k<Jz3D->dims_[2] ; k++ ) {
                        ( *Jz3D )( i, j, k ) += ( *Jz3D_s )( i, j, k );
                    }
        }
        if( rho_s[ispec] ) {
            Field3D *rho3D_s  = static_cast<Field3D *>( rho_s[ispec] );
            for( unsigned int i=0 ; i<rho3D->dims_[0] ; i++ )
                for( unsigned int j=0 ; j<rho3D->dims_[1] ; j++ )
                    for( unsigned int k=0 ; k<rho3D->dims_[2] ; k++ ) {
                        ( *rho3D )( i, j, k ) += ( *rho3D_s )( i, j, k );
                    }
        }
        
    }//END loop on species ispec
//END computeTotalRhoJ
}

// ---------------------------------------------------------------------------------------------------------------------
// Compute the total susceptibility from species susceptibility
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computeTotalEnvChi()
{
    // static cast of the total susceptibility
    Field3D *Env_Chi3D   = static_cast<Field3D *>( Env_Chi_ );
    
    // -----------------------------------
    // Species susceptibility
    // -----------------------------------
    for( unsigned int ispec=0; ispec<n_species; ispec++ ) {
        if( Env_Chi_s[ispec] ) {
            Field3D *Env_Chi3D_s  = static_cast<Field3D *>( Env_Chi_s[ispec] );
            for( unsigned int i=0 ; i<nx_p ; i++ )
                for( unsigned int j=0 ; j<ny_p ; j++ )
                    for( unsigned int k=0 ; k<nz_p ; k++ ) {
                        ( *Env_Chi3D )( i, j, k ) += ( *Env_Chi3D_s )( i, j, k );
                    }
        }
    }//END loop on species ispec
} //END computeTotalEnvChi

// ---------------------------------------------------------------------------------------------------------------------
// Compute electromagnetic energy flows vectors on the border of the simulation box
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagn3D::computePoynting()
{

    Field3D *Ex3D     = static_cast<Field3D *>( Ex_ );
    Field3D *Ey3D     = static_cast<Field3D *>( Ey_ );
    Field3D *Ez3D     = static_cast<Field3D *>( Ez_ );
    Field3D *Bx3D_m   = static_cast<Field3D *>( Bx_m );
    Field3D *By3D_m   = static_cast<Field3D *>( By_m );
    Field3D *Bz3D_m   = static_cast<Field3D *>( Bz_m );
    
    if( isXmin ) {
        unsigned int iEy=istart[0][Ey3D->isDual( 0 )];
        unsigned int iBz=istart[0][Bz3D_m->isDual( 0 )];
        unsigned int iEz=istart[0][Ez3D->isDual( 0 )];
        unsigned int iBy=istart[0][By3D_m->isDual( 0 )];
        
        unsigned int jEy=istart[1][Ey3D->isDual( 1 )];
        unsigned int jBz=istart[1][Bz3D_m->isDual( 1 )];
        unsigned int jEz=istart[1][Ez3D->isDual( 1 )];
        unsigned int jBy=istart[1][By3D_m->isDual( 1 )];
        
        unsigned int kEy=istart[2][Ey3D->isDual( 2 )];
        unsigned int kBz=istart[2][Bz3D_m->isDual( 2 )];
        unsigned int kEz=istart[2][Ez3D->isDual( 2 )];
//         unsigned int kBy=istart[2][By3D_m->isDual(2)];

        poynting_inst[0][0] = 0.;
        for( unsigned int j=0; j<=bufsize[1][Ez3D->isDual( 1 )]; j++ ) {
            for( unsigned int k=0; k<=bufsize[2][Ey3D->isDual( 2 )]; k++ ) {
            
                double Ey__ = 0.5 *( ( *Ey3D )( iEy, jEy+j,   kEy+k )   + ( *Ey3D )( iEy,   jEy+j+1, kEy+k ) );
                double Bz__ = 0.25*( ( *Bz3D_m )( iBz, jBz+j,   kBz+k )   + ( *Bz3D_m )( iBz+1, jBz+j,   kBz+k )
                                     +( *Bz3D_m )( iBz, jBz+j+1, kBz+k )   + ( *Bz3D_m )( iBz+1, jBz+j+1, kBz+k ) );
                double Ez__ = 0.5 *( ( *Ez3D )( iEz, jEz+j, kEz+k )   + ( *Ez3D )( iEz,   jEz+j,   kEz+k+1 ) );
                double By__ = 0.25*( ( *By3D_m )( iBy, jBy+j,   kEz+k )   + ( *By3D_m )( iBy+1, jBy+j,   kEz+k )
                                     +( *By3D_m )( iBy, jBy+j,   kEz+k+1 ) + ( *By3D_m )( iBy+1, jBy+j,   kEz+k+1 ) );
                                     
                poynting_inst[0][0] += Ey__*Bz__ - Ez__*By__;
            }
        }
        poynting_inst[0][0] *= dy*dz*timestep;
        poynting[0][0]+= poynting_inst[0][0];
    }//if Xmin
    
    
    if( isXmax ) {
        unsigned int offset = bufsize[0][Ey3D  ->isDual( 0 )];
        
        unsigned int iEy=istart[0][Ey3D  ->isDual( 0 )] + offset;
        unsigned int iBz=istart[0][Bz3D_m->isDual( 0 )] + offset;
        unsigned int iEz=istart[0][Ez3D  ->isDual( 0 )] + offset;
        unsigned int iBy=istart[0][By3D_m->isDual( 0 )] + offset;
        
        unsigned int jEy=istart[1][Ey3D  ->isDual( 1 )];
        unsigned int jBz=istart[1][Bz3D_m->isDual( 1 )];
        unsigned int jEz=istart[1][Ez3D  ->isDual( 1 )];
        unsigned int jBy=istart[1][By3D_m->isDual( 1 )];
        
        unsigned int kEy=istart[2][Ey3D  ->isDual( 2 )];
        unsigned int kBz=istart[2][Bz3D_m->isDual( 2 )];
        unsigned int kEz=istart[2][Ez3D  ->isDual( 2 )];
        //unsigned int kBy=istart[2][By3D_m->isDual(2)];
        
        poynting_inst[1][0] = 0.;
        for( unsigned int j=0; j<=bufsize[1][Ez3D->isDual( 1 )]; j++ ) {
            for( unsigned int k=0; k<=bufsize[2][Ey3D->isDual( 2 )]; k++ ) {
            
                double Ey__ = 0.5 *( ( *Ey3D )( iEy, jEy+j,   kEy+k )   + ( *Ey3D )( iEy,   jEy+j+1, kEy+k ) );
                double Bz__ = 0.25*( ( *Bz3D_m )( iBz, jBz+j,   kBz+k )   + ( *Bz3D_m )( iBz+1, jBz+j,   kBz+k )
                                     +( *Bz3D_m )( iBz, jBz+j+1, kBz+k )   + ( *Bz3D_m )( iBz+1, jBz+j+1, kBz+k ) );
                double Ez__ = 0.5 *( ( *Ez3D )( iEz, jEz+j,   kEz+k )   + ( *Ez3D )( iEz,   jEz+j,   kEz+k+1 ) );
                double By__ = 0.25*( ( *By3D_m )( iBy, jBy+j,   kEz+k )   + ( *By3D_m )( iBy+1, jBy+j,   kEz+k )
                                     +( *By3D_m )( iBy, jBy+j,   kEz+k+1 ) + ( *By3D_m )( iBy+1, jBy+j,   kEz+k+1 ) );
                                     
                poynting_inst[1][0] += Ey__*Bz__ - Ez__*By__;
            }
        }
        poynting_inst[1][0] *= dy*dz*timestep;
        poynting[1][0] -= poynting_inst[1][0];
    }//if Xmax
    
    if( isYmin ) {
    
        unsigned int iEz=istart[0][Ez_->isDual( 0 )];
        unsigned int iBx=istart[0][Bx_m->isDual( 0 )];
        unsigned int iEx=istart[0][Ex_->isDual( 0 )];
        unsigned int iBz=istart[0][Bz_m->isDual( 0 )];
        
        unsigned int jEz=istart[1][Ez_->isDual( 1 )];
        unsigned int jBx=istart[1][Bx_m->isDual( 1 )];
        unsigned int jEx=istart[1][Ex_->isDual( 1 )];
        unsigned int jBz=istart[1][Bz_m->isDual( 1 )];
        
        unsigned int kEz=istart[2][Ez_->isDual( 2 )];
        unsigned int kBx=istart[2][Bx_m->isDual( 2 )];
        unsigned int kEx=istart[2][Ex_->isDual( 2 )];
//         unsigned int kBz=istart[2][Bz_m->isDual(2)];

        poynting_inst[0][1] = 0.;
        for( unsigned int i=0; i<=bufsize[0][Ez3D->isDual( 0 )]; i++ ) {
            for( unsigned int k=0; k<=bufsize[2][Ex3D->isDual( 2 )]; k++ ) {
                double Ez__ = 0.5 *( ( *Ez3D )( iEz+i, jEz,   kEz+k )   + ( *Ez3D )( iEz+i,   jEz,   kEz+k+1 ) );
                double Bx__ = 0.25*( ( *Bx3D_m )( iBx+i, jBx,   kBx+k )   + ( *Bx3D_m )( iBx+i,   jBx+1, kBx+k )
                                     +( *Bx3D_m )( iBx+i, jBx,   kBx+k+1 ) + ( *Bx3D_m )( iBx+i,   jBx+1, kBx+k+1 ) );
                double Ex__ = 0.5 *( ( *Ex3D )( iEx+i, jEx,   kEx+k )   + ( *Ex3D )( iEx+i+1, jEx,   kEx+k ) );
                double Bz__ = 0.25*( ( *Bz3D_m )( iBz+i, jBz,   kEx+k )   + ( *Bz3D_m )( iBz+i+1, jBz,   kEx+k )
                                     +( *Bz3D_m )( iBz+i, jBz+1, kEx+k )   + ( *Bz3D_m )( iBz+i+1, jBz+1, kEx+k ) );
                                     
                poynting_inst[0][1] += Ez__*Bx__ - Ex__*Bz__;
            }
        }
        poynting_inst[0][1] *= dx*dz*timestep;
        poynting[0][1] += poynting_inst[0][1];
    }// if Ymin
    
    if( isYmax ) {
    
        unsigned int iEz=istart[0][Ez3D  ->isDual( 0 )];
        unsigned int iBx=istart[0][Bx3D_m->isDual( 0 )];
        unsigned int iEx=istart[0][Ex3D  ->isDual( 0 )];
        unsigned int iBz=istart[0][Bz3D_m->isDual( 0 )];
        
        unsigned int offset = bufsize[1][Ez3D->isDual( 1 )];
        
        unsigned int jEz=istart[1][Ez3D  ->isDual( 1 )] + offset;
        unsigned int jBx=istart[1][Bx3D_m->isDual( 1 )] + offset;
        unsigned int jEx=istart[1][Ex3D  ->isDual( 1 )] + offset;
        unsigned int jBz=istart[1][Bz3D_m->isDual( 1 )] + offset;
        
        unsigned int kEz=istart[2][Ez_  ->isDual( 2 )];
        unsigned int kBx=istart[2][Bx_m->isDual( 2 )];
        unsigned int kEx=istart[2][Ex_  ->isDual( 2 )];
        //unsigned int kBz=istart[2][Bz_m->isDual(2)];
        
        poynting_inst[1][1] = 0.;
        for( unsigned int i=0; i<=bufsize[0][Ez3D->isDual( 0 )]; i++ ) {
            for( unsigned int k=0; k<=bufsize[2][Ex3D->isDual( 2 )]; k++ ) {
                double Ez__ = 0.5 *( ( *Ez3D )( iEz+i, jEz,   kEz+k )   + ( *Ez3D )( iEz+i,   jEz,   kEz+k+1 ) );
                double Bx__ = 0.25*( ( *Bx3D_m )( iBx+i, jBx,   kBx+k )   + ( *Bx3D_m )( iBx+i,   jBx+1, kBx+k )
                                     +( *Bx3D_m )( iBx+i, jBx,   kBx+k+1 ) + ( *Bx3D_m )( iBx+i,   jBx+1, kBx+k+1 ) );
                double Ex__ = 0.5 *( ( *Ex3D )( iEx+i, jEx,   kEx+k )   + ( *Ex3D )( iEx+i+1, jEx,   kEx+k ) );
                double Bz__ = 0.25*( ( *Bz3D_m )( iBz+i, jBz,   kEx+k )   + ( *Bz3D_m )( iBz+i+1, jBz,   kEx+k )
                                     +( *Bz3D_m )( iBz+i, jBz+1, kEx+k )   + ( *Bz3D_m )( iBz+i+1, jBz+1, kEx+k ) );
                                     
                poynting_inst[1][1] += Ez__*Bx__ - Ex__*Bz__;
            }
        }
        poynting_inst[1][1] *= dx*dz*timestep;
        poynting[1][1] -= poynting_inst[1][1];
    }//if Ymax
    
    if( isZmin ) {
    
        unsigned int iEx=istart[0][Ex_->isDual( 0 )];
        unsigned int iBy=istart[0][By_m->isDual( 0 )];
        unsigned int iEy=istart[0][Ey_->isDual( 0 )];
        unsigned int iBx=istart[0][Bx_m->isDual( 0 )];
        
        unsigned int jEx=istart[1][Ex_->isDual( 1 )];
        unsigned int jBy=istart[1][By_m->isDual( 1 )];
        unsigned int jEy=istart[1][Ey_->isDual( 1 )];
        unsigned int jBx=istart[1][Bx_m->isDual( 1 )];
        
        unsigned int kEx=istart[2][Ex_->isDual( 2 )];
        unsigned int kBy=istart[2][By_m->isDual( 2 )];
        unsigned int kEy=istart[2][Ey_->isDual( 2 )];
//         unsigned int kBx=istart[2][Bx_m->isDual(2)];

        poynting_inst[0][2] = 0.;
        for( unsigned int i=0; i<=bufsize[0][Ez3D->isDual( 0 )]; i++ ) {
            for( unsigned int j=0; j<=bufsize[1][Ex3D->isDual( 1 )]; j++ ) {
            
                double Ex__ = 0.5 *( ( *Ex3D )( iEx+i, jEx+j, kEx )   + ( *Ex3D )( iEx+i+1, jEx+j,   kEx ) );
                double By__ = 0.25*( ( *By3D_m )( iBy+i, jBy+j, kBy )   + ( *By3D_m )( iBy+i+1, jBy+j,   kBy )
                                     +( *By3D_m )( iBy+i, jBy+j, kBy+1 ) + ( *By3D_m )( iBy+i+1, jBy+j,   kBy+1 ) );
                double Ey__ = 0.5 *( ( *Ey3D )( iEy+i, jEy+j, kEy )   + ( *Ey3D )( iEy+i,   jEy+j+1, kEy ) );
                double Bx__ = 0.25*( ( *Bx3D_m )( iBx+i, jBx+j, kEx )   + ( *Bx3D_m )( iBx+i,   jBx+j+1, kEx )
                                     +( *Bx3D_m )( iBx+i, jBx+j, kEx+1 ) + ( *Bx3D_m )( iBx+i,   jBx+j+1, kEx+1 ) );
                                     
                poynting_inst[0][2] += Ex__*By__ - Ey__*Bx__;
            }
        }
        poynting_inst[0][2] *= dx*dy*timestep;
        poynting[0][2] += poynting_inst[0][2];
    }
    
    if( isZmax ) {
    
        unsigned int iEx=istart[0][Ex_ ->isDual( 0 )];
        unsigned int iBy=istart[0][By_m->isDual( 0 )];
        unsigned int iEy=istart[0][Ey_ ->isDual( 0 )];
        unsigned int iBx=istart[0][Bx_m->isDual( 0 )];
        
        unsigned int jEx=istart[1][Ex_ ->isDual( 1 )];
        unsigned int jBy=istart[1][By_m->isDual( 1 )];
        unsigned int jEy=istart[1][Ey_ ->isDual( 1 )];
        unsigned int jBx=istart[1][Bx_m->isDual( 1 )];
        
        unsigned int offset = bufsize[2][Ex_->isDual( 2 )];
        
        unsigned int kEx=istart[2][Ex_ ->isDual( 2 )] + offset;
        unsigned int kBy=istart[2][By_m->isDual( 2 )] + offset;
        unsigned int kEy=istart[2][Ey_ ->isDual( 2 )] + offset;
//        unsigned int kBx=istart[2][Bx_m->isDual(2)] + offset;

        poynting_inst[1][2] = 0.;
        for( unsigned int i=0; i<=bufsize[0][Ez3D->isDual( 0 )]; i++ ) {
            for( unsigned int j=0; j<=bufsize[1][Ex3D->isDual( 1 )]; j++ ) {
            
                double Ex__ = 0.5 *( ( *Ex3D )( iEx+i, jEx+j, kEx )   + ( *Ex3D )( iEx+i+1, jEx+j,   kEx ) );
                double By__ = 0.25*( ( *By3D_m )( iBy+i, jBy+j, kBy )   + ( *By3D_m )( iBy+i+1, jBy+j,   kBy )
                                     +( *By3D_m )( iBy+i, jBy+j, kBy+1 ) + ( *By3D_m )( iBy+i+1, jBy+j,   kBy+1 ) );
                double Ey__ = 0.5 *( ( *Ey3D )( iEy+i, jEy+j, kEy )   + ( *Ey3D )( iEy+i,   jEy+j+1, kEy ) );
                double Bx__ = 0.25*( ( *Bx3D_m )( iBx+i, jBx+j, kEx )   + ( *Bx3D_m )( iBx+i,   jBx+j+1, kEx )
                                     +( *Bx3D_m )( iBx+i, jBx+j, kEx+1 ) + ( *Bx3D_m )( iBx+i,   jBx+j+1, kEx+1 ) );
                                     
                poynting_inst[1][2] += Ex__*By__ - Ey__*Bx__;
            }
        }
        poynting_inst[1][2] *= dx*dy*timestep;
        poynting[1][2] -= poynting_inst[1][2];
    }
    
}

void ElectroMagn3D::applyExternalField( Field *my_field,  Profile *profile, Patch *patch )
{

    Field3D *field3D = static_cast<Field3D *>( my_field );
    
    vector<bool> dual(3, false);
    if ( ( field3D->name.substr(0,2) == "Jx" ) || ( field3D->name.substr(0,2) == "Ex" ) )
        dual[0] = true;
    if ( ( field3D->name.substr(0,2) == "Jy" ) || ( field3D->name.substr(0,2) == "Ey" ) )
        dual[1] = true;
    if ( ( field3D->name.substr(0,2) == "Jz" ) || ( field3D->name.substr(0,2) == "Ez" ) )
        dual[2] = true;

    if ( field3D->name.substr(0,2) == "Bx" ) {
        dual[1] = true;
        dual[2] = true;
    }
    if ( field3D->name.substr(0,2) == "By" ) {
        dual[0] = true;
        dual[2] = true;
    }
    if ( field3D->name.substr(0,2) == "Bz" ) {
        dual[0] = true;
        dual[1] = true;
    }

    vector<double> pos( 3 );
    pos[0]      = dx*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( dual[0]?-0.5:0. ) );
    double pos1 = dy*( ( double )( patch->getCellStartingGlobalIndex( 1 ) )+( dual[1]?-0.5:0. ) );
    double pos2 = dz*( ( double )( patch->getCellStartingGlobalIndex( 2 ) )+( dual[2]?-0.5:0. ) );
    int N0 = ( int )field3D->dims()[0];
    int N1 = ( int )field3D->dims()[1];
    int N2 = ( int )field3D->dims()[2];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    // Create the x,y,z maps where profiles will be evaluated
    vector<Field *> xyz( 3 );
    vector<unsigned int> n_space_to_create( 3 );
    n_space_to_create[0] = N0;
    n_space_to_create[1] = N1;
    n_space_to_create[2] = N2;
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        xyz[idim] = new Field3D( n_space_to_create );
    }
    
    for( int i=0 ; i<N0 ; i++ ) {
        pos[1] = pos1;
        for( int j=0 ; j<N1 ; j++ ) {
            pos[2] = pos2;
            for( int k=0 ; k<N2 ; k++ ) {
                //(*field3D)(i,j,k) += profile->valueAt(pos);
                for( unsigned int idim=0 ; idim<3 ; idim++ ) {
                    ( *xyz[idim] )( i, j, k ) = pos[idim];
                }
                pos[2] += dz;
            }
            pos[1] += dy;
        }
        pos[0] += dx;
    }
    
    profile->valuesAt( xyz, *field3D );
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        delete xyz[idim];
    }
    
    
}

void ElectroMagn3D::applyPrescribedField( Field *my_field,  Profile *profile, Patch *patch, double time )
{

    Field3D *field3D = static_cast<Field3D *>( my_field );
    
    vector<double> pos( 3 );
    pos[0]      = dx*( ( double )( patch->getCellStartingGlobalIndex( 0 ) )+( field3D->isDual( 0 )?-0.5:0. ) );
    double pos1 = dy*( ( double )( patch->getCellStartingGlobalIndex( 1 ) )+( field3D->isDual( 1 )?-0.5:0. ) );
    double pos2 = dz*( ( double )( patch->getCellStartingGlobalIndex( 2 ) )+( field3D->isDual( 2 )?-0.5:0. ) );
    int N0 = ( int )field3D->dims()[0];
    int N1 = ( int )field3D->dims()[1];
    int N2 = ( int )field3D->dims()[2];
    
    // UNSIGNED INT LEADS TO PB IN PERIODIC BCs
    // Create the x,y,z maps where profiles will be evaluated
    vector<Field *> xyz( 3 );
    vector<unsigned int> n_space_to_create( 3 );
    n_space_to_create[0] = N0;
    n_space_to_create[1] = N1;
    n_space_to_create[2] = N2;
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        xyz[idim] = new Field3D( n_space_to_create );
    }
    
    for( int i=0 ; i<N0 ; i++ ) {
        pos[1] = pos1;
        for( int j=0 ; j<N1 ; j++ ) {
            pos[2] = pos2;
            for( int k=0 ; k<N2 ; k++ ) {
                //(*field3D)(i,j,k) += profile->valueAt(pos);
                for( unsigned int idim=0 ; idim<3 ; idim++ ) {
                    ( *xyz[idim] )( i, j, k ) = pos[idim];
                }
                pos[2] += dz;
            }
            pos[1] += dy;
        }
        pos[0] += dx;
    }
    
    profile->valuesAtTime( xyz, time, *field3D );
    
    for( unsigned int idim=0 ; idim<3 ; idim++ ) {
        delete xyz[idim];
    }

}



void ElectroMagn3D::initAntennas( Patch *patch, Params& params )
{

    // Filling the space profiles of antennas
    for( unsigned int i=0; i<antennas.size(); i++ ) {
        if      (antennas[i].fieldName == "Jx")
            antennas[i].field = FieldFactory::create(dimPrim, 0, false, "Jx", params);
        else if (antennas[i].fieldName == "Jy")
            antennas[i].field = FieldFactory::create(dimPrim, 1, false, "Jy", params);
        else if (antennas[i].fieldName == "Jz")
            antennas[i].field = FieldFactory::create(dimPrim, 2, false, "Jz", params);
        else {
            ERROR("Antenna cannot be applied to field "<<antennas[i].fieldName);
        }
        
        if( antennas[i].field ) {
            applyExternalField( antennas[i].field, antennas[i].space_profile, patch );
        }
    }
    
}

