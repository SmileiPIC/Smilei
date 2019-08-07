
#include "MA_Solver3D_norm.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MA_Solver3D_norm::MA_Solver3D_norm( Params &params )
    : Solver3D( params )
{
}

MA_Solver3D_norm::~MA_Solver3D_norm()
{
}

void MA_Solver3D_norm::operator()( ElectroMagn *fields )
{

    // Static-cast of the fields
    Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( fields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
    Field3D *Jx3D = static_cast<Field3D *>( fields->Jx_ );
    Field3D *Jy3D = static_cast<Field3D *>( fields->Jy_ );
    Field3D *Jz3D = static_cast<Field3D *>( fields->Jz_ );
    
    // Electric field Ex^(d,p,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Ex3D )( i, j, k ) += -dt*( *Jx3D )( i, j, k )
                                        +                 dt_ov_dy * ( ( *Bz3D )( i, j+1, k ) - ( *Bz3D )( i, j, k ) )
                                        -                 dt_ov_dz * ( ( *By3D )( i, j, k+1 ) - ( *By3D )( i, j, k ) );
            }
        }
    }
    
    // Electric field Ey^(p,d,p)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Ey3D )( i, j, k ) += -dt*( *Jy3D )( i, j, k )
                                        -                  dt_ov_dx * ( ( *Bz3D )( i+1, j, k ) - ( *Bz3D )( i, j, k ) )
                                        +                  dt_ov_dz * ( ( *Bx3D )( i, j, k+1 ) - ( *Bx3D )( i, j, k ) );
            }
        }
    }
    
    // Electric field Ez^(p,p,d)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_d ; k++ ) {
                ( *Ez3D )( i, j, k ) += -dt*( *Jz3D )( i, j, k )
                                        +                  dt_ov_dx * ( ( *By3D )( i+1, j, k ) - ( *By3D )( i, j, k ) )
                                        -                  dt_ov_dy * ( ( *Bx3D )( i, j+1, k ) - ( *Bx3D )( i, j, k ) );
            }
        }
    }
    
}

