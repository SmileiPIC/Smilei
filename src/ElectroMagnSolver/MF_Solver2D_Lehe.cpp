
#include "MF_Solver2D_Lehe.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include <algorithm>

MF_Solver2D_Lehe::MF_Solver2D_Lehe( Params &params )
    : Solver2D( params )
{
    //ERROR("Under development, not yet working");
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    
    beta_yx = 1./8.;
    beta_xy = pow( dx/dy, 2 )/8.;
    delta_x = ( 1./4. )*( 1.-pow( sin( M_PI*dt_ov_dx/2. )/dt_ov_dx, 2 ) );
    
    alpha_y =  1.-2.*beta_yx;
    alpha_x =  1.-2.*beta_xy-3.*delta_x;
}

MF_Solver2D_Lehe::~MF_Solver2D_Lehe()
{
}

void MF_Solver2D_Lehe::operator()( ElectroMagn *fields )
{
    // Static-cast of the fields
    Field2D *Ex2D = static_cast<Field2D *>( fields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( fields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    
    
    // Magnetic field Bx^(p,d)
    for( unsigned int i=1 ; i<nx_p-1;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bx2D )( i, j ) -= dt_ov_dy * ( alpha_y * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i, j-1 ) )
                                              + beta_yx * ( ( *Ez2D )( i+1, j )-( *Ez2D )( i+1, j-1 ) + ( *Ez2D )( i-1, j )-( *Ez2D )( i-1, j-1 ) ) );
        }
    }
    
    // Magnetic field By^(d,p)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_p-1 ; j++ ) {
            ( *By2D )( i, j ) += dt_ov_dx * ( alpha_x * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i-1, j ) )
                                              + beta_xy * ( ( *Ez2D )( i, j+1 )-( *Ez2D )( i-1, j+1 ) + ( *Ez2D )( i, j-1 )-( *Ez2D )( i-1, j-1 ) )
                                              + delta_x * ( ( *Ez2D )( i+1, j ) - ( *Ez2D )( i-2, j ) ) );
        }
    }
    
    
    // Magnetic field Bz^(d,d)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bz2D )( i, j ) += dt_ov_dy * ( alpha_y * ( ( *Ex2D )( i, j )-( *Ex2D )( i, j-1 ) )
                                              + beta_yx * ( ( *Ex2D )( i+1, j )-( *Ex2D )( i+1, j-1 ) + ( *Ex2D )( i-1, j )-( *Ex2D )( i-1, j-1 ) ) )
                                 - dt_ov_dx * ( alpha_x * ( ( *Ey2D )( i, j )-( *Ey2D )( i-1, j ) )
                                                + beta_xy * ( ( *Ey2D )( i, j+1 )-( *Ey2D )( i-1, j+1 ) + ( *Ey2D )( i, j-1 )-( *Ey2D )( i-1, j-1 ) )
                                                + delta_x * ( ( *Ey2D )( i+1, j )-( *Ey2D )( i-2, j ) ) );
        }
    }
    
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    for( unsigned int j=0 ; j<ny_p ; j++ ) {
        ( *By2D )( 1, j ) += dt_ov_dx * ( ( *Ez2D )( 1, j ) - ( *Ez2D )( 0, j ) );
    }
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    for( unsigned int j=0 ; j<ny_p ; j++ ) {
        ( *By2D )( nx_d-2, j ) += dt_ov_dx * ( ( *Ez2D )( nx_d-2, j ) - ( *Ez2D )( nx_d-3, j ) );
    }
    
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
        ( *Bz2D )( 1, j ) += dt_ov_dx * ( ( *Ey2D )( 0, j ) - ( *Ey2D )( 1, j ) )
                             +               dt_ov_dy * ( ( *Ex2D )( 1, j ) - ( *Ex2D )( 1, j-1 ) );
    }
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {
        ( *Bz2D )( nx_d-2, j ) += dt_ov_dx * ( ( *Ey2D )( nx_d-3, j ) - ( *Ey2D )( nx_d-2, j ) )
                                  +                    dt_ov_dy * ( ( *Ex2D )( nx_d-2, j ) - ( *Ex2D )( nx_d-2, j-1 ) );
    }
    
//}// end parallel
}//END solveMaxwellFaraday



