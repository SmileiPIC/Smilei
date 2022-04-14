#include "MF_Solver3D_M4.h"

#include "ElectroMagn.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

#include <algorithm>

MF_Solver3D_M4::MF_Solver3D_M4( Params &params )
    : Solver3D( params )
{

    beta_xy = dt_ov_dy * dt_ov_dy / 12.0;
    beta_yx = dt_ov_dx * dt_ov_dx / 12.0;
    beta_xz = dt_ov_dz * dt_ov_dz / 12.0;
    beta_zx = beta_yx;
    beta_yz = beta_xz;
    beta_zy = beta_xy;
    delta_x = beta_yx - 1.0 / 12.0;
    delta_y = beta_xy - 1.0 / 12.0;
    delta_z = beta_xz - 1.0 / 12.0;
    alpha_x = 1.0 - 2.0 * ( beta_xy+beta_xz ) - 3.0 * delta_x;
    alpha_y = 1.0 - 2.0 * ( beta_yx+beta_yz ) - 3.0 * delta_y;
    alpha_z = 1.0 - 2.0 * ( beta_zx+beta_zy ) - 3.0 * delta_z;

    Ax  = alpha_x*dt_ov_dx;
    Ay  = alpha_y*dt_ov_dy;
    Az  = alpha_z*dt_ov_dz;
    Bxy = beta_xy*dt_ov_dx;
    Byx = beta_yx*dt_ov_dy;
    Bxz = beta_xz*dt_ov_dx;
    Bzx = beta_zx*dt_ov_dz;
    Byz = beta_yz*dt_ov_dy;
    Bzy = beta_zy*dt_ov_dz;
    Dx  = delta_x*dt_ov_dx;
    Dy  = delta_y*dt_ov_dy;
    Dz  = delta_z*dt_ov_dz;

}

MF_Solver3D_M4::~MF_Solver3D_M4()
{
}

void MF_Solver3D_M4::operator()( ElectroMagn *fields )
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
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                ( *Bx3D )( i, j, k ) += Az * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i, j, k-1 ) )
                                     + Bzx * ( ( *Ey3D )( i+1, j, k )-( *Ey3D )( i+1, j, k-1 ) + ( *Ey3D )( i-1, j, k )-( *Ey3D )( i-1, j, k-1 ) )
                                     + Bzy * ( ( *Ey3D )( i, j+1, k )-( *Ey3D )( i, j+1, k-1 ) + ( *Ey3D )( i, j-1, k )-( *Ey3D )( i, j-1, k-1 ) )
                                     +  Dz * ( ( *Ey3D )( i, j, k+1 ) - ( *Ey3D )( i, j, k-2) )
                                     -  Ay * ( ( *Ez3D )( i,  j, k )  - ( *Ez3D )( i,  j-1, k ) )
                                     - Byx * ( ( *Ez3D )( i+1, j, k )  - ( *Ez3D )( i+1, j-1, k ) + ( *Ez3D )( i-1, j, k )-( *Ez3D )( i-1, j-1, k ) )
                                     - Byz * ( ( *Ez3D )( i, j, k+1 )  - ( *Ez3D )( i, j-1, k+1 ) + ( *Ez3D )( i, j, k-1 )-( *Ez3D )( i, j-1, k-1 ) )
                                     -  Dy * ( ( *Ez3D )( i, j+1, k )-( *Ez3D )( i, j-2, k ) );

            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
            for( unsigned int k=2 ; k<nz_d-2 ; k++ ) {
                ( *By3D )( i, j, k ) += Ax * ( ( *Ez3D )( i,  j, k ) - ( *Ez3D )( i-1, j, k ) )
                                     + Bxy * ( ( *Ez3D )( i,  j+1, k ) - ( *Ez3D )( i-1, j+1, k ) + ( *Ez3D )( i, j-1, k )-( *Ez3D )( i-1, j-1, k ) )
                                     + Bxz * ( ( *Ez3D )( i, j, k+1 ) - ( *Ez3D )( i-1, j, k+1 ) + ( *Ez3D )( i, j, k-1 )-( *Ez3D )( i-1, j, k-1 ) )
                                     +  Dx * ( ( *Ez3D )( i+1, j, k ) - ( *Ez3D )( i-2, j, k ) ) 
                                     -  Az * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j, k-1 ) )
                                     - Bzx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j, k-1 ) + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j, k-1 ) )
                                     - Bzy * ( ( *Ex3D )( i, j+1, k )-( *Ex3D )( i, j+1, k-1 ) + ( *Ex3D )( i, j-1, k )-( *Ex3D )( i, j-1, k-1 ) )
                                     -  Dz * ( ( *Ex3D )( i, j, k+1 ) - ( *Ex3D )( i, j, k-2) ) ;
            }
        }
    } 
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
            for( unsigned int k=1 ; k<nz_p-1 ; k++ ) {
                ( *Bz3D )( i, j, k ) += Ay * ( ( *Ex3D )( i, j, k )-( *Ex3D )( i, j-1, k ) )
                                     + Byx * ( ( *Ex3D )( i+1, j, k )-( *Ex3D )( i+1, j-1, k ) + ( *Ex3D )( i-1, j, k )-( *Ex3D )( i-1, j-1, k )) 
                                     + Byz * ( ( *Ex3D )( i, j, k+1 )-( *Ex3D )( i, j-1, k+1 ) + ( *Ex3D )( i, j, k-1 )-( *Ex3D )( i, j-1, k-1 ))
                                     + Dy  * ( ( *Ex3D )( i, j+1, k )-( *Ex3D )( i, j-2, k ) )
                                     -  Ax * ( ( *Ey3D )( i, j, k )-( *Ey3D )( i-1, j, k ) )
                                     - Bxy * ( ( *Ey3D )( i, j+1, k )-( *Ey3D )( i-1, j+1, k ) + ( *Ey3D )( i, j-1, k )-( *Ey3D )( i-1, j-1, k ))
                                     - Bxz * ( ( *Ey3D )( i, j, k+1 )-( *Ey3D )( i-1, j, k+1 ) + ( *Ey3D )( i, j, k-1 )-( *Ey3D )( i-1, j, k-1 ))
                                     - Dx  * ( ( *Ey3D )( i+1, j, k )-( *Ey3D )( i-2, j, k ) ) ;
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



