
#include "MF_Solver3D_Lehe.h"

#include "ElectroMagn.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

#include <algorithm>

MF_Solver3D_Lehe::MF_Solver3D_Lehe( Params &params )
    : Solver3D( params )
{
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    
    beta_yx = 1./8.; // = beta_zx as well but we define and use only 1 variable
    beta_xy = pow( dx/dy, 2 )/8.;
    beta_xz = pow( dx/dz, 2 )/8.;
    delta_x = ( 1./4. )*( 1.-pow( sin( M_PI*dt_ov_dx/2. )/dt_ov_dx, 2 ) );
    
    alpha_y =  1. - 2.*beta_yx; // = alpha_z as well but we define and use only 1 variable
    alpha_x =  1. - 2.*beta_xy - 2.*beta_xz - 3.*delta_x ;
    
    
}

MF_Solver3D_Lehe::~MF_Solver3D_Lehe()
{
}

void MF_Solver3D_Lehe::operator()( ElectroMagn *fields )
{
    // Static-cast of the fields
    Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( fields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
    ElectroMagn3D *EM3D = static_cast<ElectroMagn3D *>( fields );
    
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=1 ; i<nx_p-1;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *Bx3D )( i, j, k ) += -dt_ov_dy * ( alpha_y * ( ( *Ez3D )( i,  j, k )  - ( *Ez3D )( i,  j-1, k ) )
                                                      + beta_yx * ( ( *Ez3D )( i+1, j, k )  - ( *Ez3D )( i+1, j-1, k ) + ( *Ez3D )( i-1, j, k )-( *Ez3D )( i-1, j-1, k ) )
                                                    )
                                        +dt_ov_dz * ( /*alpha_z*/ alpha_y * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i, j, k-1 ) )
                                                +/*beta_zx*/ beta_yx * ( ( *Ey3D )( i+1, j, k )-( *Ey3D )( i+1, j, k-1 ) + ( *Ey3D )( i-1, j, k )-( *Ey3D )( i-1, j, k-1 ) )
                                                    );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( i, j, k ) += dt_ov_dx * ( alpha_x * ( ( *Ez3D )( i,  j, k ) - ( *Ez3D )( i-1, j, k ) )
                                                     + beta_xy * ( ( *Ez3D )( i,  j+1, k ) - ( *Ez3D )( i-1, j+1, k ) + ( *Ez3D )( i, j-1, k )-( *Ez3D )( i-1, j-1, k ) )
                                                     + beta_xz * ( ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i-1, j, k+1 ) + ( *Ez3D )( i, j, k-1 )-( *Ez3D )( i-1, j, k-1 ) )
                                                     + delta_x * ( ( *Ez3D )( i+1, j, k ) - ( *Ez3D )( i-2, j, k ) ) )
                                        -dt_ov_dz * ( /*alpha_z*/ alpha_y * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j, k-1 ) )
                                                +/*beta_zx*/ beta_yx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j, k-1 ) + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j, k-1 ) )
                                                    );
            }
        }
    }
    
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_p-1 ; k++ ) {
                ( *Bz3D )( i, j, k ) += dt_ov_dy * ( alpha_y * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j-1, k ) )
                                                     + beta_yx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j-1, k )   + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j-1, k ) ) )
                                        - dt_ov_dx * ( alpha_x * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i-1, j, k ) )
                                                       + beta_xy * ( ( *Ey3D )( i, j+1, k )-( *Ey3D )( i-1, j+1, k )   + ( *Ey3D )( i, j-1, k )-( *Ey3D )( i-1, j-1, k ) )
                                                       + beta_xz * ( ( *Ey3D )( i, j, k+1 )-( *Ey3D )( i-1, j, k+1 ) + ( *Ey3D )( i, j, k-1 )-( *Ey3D )( i-1, j, k-1 ) )
                                                       + delta_x * ( ( *Ey3D )( i+1, j, k )-( *Ey3D )( i-2, j, k ) ) );
            }
        }
    }
    
    //Additional boundaries treatment on the primal direction of each B field
    
    // Magnetic field Bx^(p,d,d)
    if( EM3D->isXmin ) {
        // At Xmin
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *Bx3D )( 0, j, k ) += -dt_ov_dy * ( ( *Ez3D )( 0, j, k ) - ( *Ez3D )( 0, j-1, k ) ) + dt_ov_dz * ( ( *Ey3D )( 0, j, k ) - ( *Ey3D )( 0, j, k-1 ) );
            }
        }
        //Additional boundaries treatment for i=1 and i=nx_d-2 for By and Bz
        // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( 1, j, k ) += dt_ov_dx * ( ( *Ez3D )( 1, j, k ) - ( *Ez3D )( 0, j, k ) )
                                        -dt_ov_dz * ( ( *Ex3D )( 1, j, k ) - ( *Ex3D )( 1, j, k-1 ) );
            }
        }
        // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( 1, j, k ) += dt_ov_dx * ( ( *Ey3D )( 0, j, k ) - ( *Ey3D )( 1, j, k ) )
                                        +  dt_ov_dy * ( ( *Ex3D )( 1, j, k ) - ( *Ex3D )( 1, j-1, k ) );
            }
        }
        
    }
    if( EM3D->isXmax ) {
        // At Xmax
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *Bx3D )( nx_p-1, j, k ) += -dt_ov_dy * ( ( *Ez3D )( nx_p-1, j, k ) - ( *Ez3D )( nx_p-1, j-1, k ) ) + dt_ov_dz * ( ( *Ey3D )( nx_p-1, j, k ) - ( *Ey3D )( nx_p-1, j, k-1 ) );
            }
        }
        // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( nx_d-2, j, k ) += dt_ov_dx * ( ( *Ez3D )( nx_d-2, j, k ) - ( *Ez3D )( nx_d-3, j, k ) )
                                             -dt_ov_dz * ( ( *Ex3D )( nx_d-2, j, k ) - ( *Ex3D )( nx_d-2, j, k-1 ) );
            }
        }
        // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( nx_d-2, j, k ) += dt_ov_dx * ( ( *Ey3D )( nx_d-3, j, k ) - ( *Ey3D )( nx_d-2, j, k ) )
                                             +  dt_ov_dy * ( ( *Ex3D )( nx_d-2, j, k ) - ( *Ex3D )( nx_d-2, j-1, k ) );
            }
        }
        
    }
    
    if( EM3D->isYmin ) {
        //At Ymin
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) { //Cases i=1 and i=nx_d-2 are treated later for all required j and k.
            unsigned int j=0 ;
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( i, j, k ) += -dt_ov_dz * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j, k-1 ) ) + dt_ov_dx * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i-1, j, k ) );
            }
        }
    }
    
    if( EM3D->isYmax ) {
        //At Ymax
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) { //Cases i=1 and i=nx_d-2 are treated later for all required j and k
            unsigned int j=ny_p-1 ;
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( i, j, k ) += -dt_ov_dz * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j, k-1 ) ) + dt_ov_dx * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i-1, j, k ) );
            }
        }
    }
    
    if( EM3D->isZmin ) {
        //At Zmin
        // Magnetic field Bz^(d,d,p)
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) { //Cases i=1 and i=nx_d-2 are treated later for all required j and k
            unsigned int k=0 ;
            for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
                ( *Bz3D )( i, j, k ) += -dt_ov_dx * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i-1, j, k ) ) + dt_ov_dy * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j-1, k ) );
            }
        }
    }
    
    if( EM3D->isZmax ) {
        //At Zmax
        // Magnetic field Bz^(d,d,p)
        for( unsigned int i=2 ; i<nx_d-2 ; i++ ) { //Cases i=1 and i=nx_d-2 are treated later for all required j and k
            unsigned int k=nz_p-1 ;
            for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
                ( *Bz3D )( i, j, k ) += -dt_ov_dx * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i-1, j, k ) ) + dt_ov_dy * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j-1, k ) );
            }
        }
    }
    
    
    
    
    
    
    
    
}//END solveMaxwellFaraday



