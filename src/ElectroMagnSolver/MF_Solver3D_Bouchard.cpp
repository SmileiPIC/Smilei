#include "MF_Solver3D_Bouchard.h"

#include "ElectroMagn.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

#include <algorithm>

MF_Solver3D_Bouchard::MF_Solver3D_Bouchard( Params &params )
    : Solver3D( params )
{
    //ERROR("Under development, not yet working");
    double dt = params.timestep;
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    dz = params.cell_length[2];
    double dx_ov_dt  = dx/dt;
    double dy_ov_dt  = dy/dt;
    double dz_ov_dt  = dz/dt;
    double dt_ov_dx  = dt/dx;
    double dt_ov_dy  = dt/dy;
    double dt_ov_dz  = dt/dz;
    if( (dx!=dy)||(dx!=dz)||(dy!=dz) ) {
        ERROR( "Bouchard solver requires the same cell-length in x, y and z directions" );
    }
    if( dx_ov_dt!=2 ) {
        WARNING( "Bouchard solver requires dx/dt = 2 (Magical Timestep)" );
    }

    double delta = 0.1222*(1-pow(2.,2))/4. ;
    double beta = -0.1727*(1-0.5*pow(2.,2)-4.*delta)/4. ;
    double alpha = 1-4.*beta-3.*delta ;
    
    delta_x = delta ;
    delta_y = delta ;
    delta_z = delta ;
    beta_xy = beta ;
    beta_yx = beta ;
    beta_xz = beta ;
    beta_zx = beta ;
    beta_yz = beta ;
    beta_zy = beta ;
    alpha_x = alpha ;
    alpha_y = alpha ;
    alpha_z = alpha ;

    Ax  = alpha_x*dt/dx;
    Ay  = alpha_y*dt/dy;
    Az  = alpha_z*dt/dz;
    Bxy = beta_xy*dt/dx;
    Byx = beta_yx*dt/dy;
    Bxz = beta_xz*dt/dx;
    Bzx = beta_zx*dt/dz;
    Byz = beta_yz*dt/dy;
    Bzy = beta_zy*dt/dz;
    Dx  = delta_x*dt/dx;
    Dy  = delta_y*dt/dy;
    Dz  = delta_z*dt/dz;

}

MF_Solver3D_Bouchard::~MF_Solver3D_Bouchard()
{
}

void MF_Solver3D_Bouchard::operator()( ElectroMagn *fields )
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



