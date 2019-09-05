
#include "MF_Solver3D_Yee.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MF_Solver3D_Yee::MF_Solver3D_Yee( Params &params )
    : Solver3D( params )
{
}

MF_Solver3D_Yee::~MF_Solver3D_Yee()
{
}

void MF_Solver3D_Yee::operator()( ElectroMagn *fields )
{
    // Static-cast of the fields
    Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( fields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *Bx3D )( i, j, k ) += -dt_ov_dy * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i, j-1, k ) ) + dt_ov_dz * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i, j, k-1 ) );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( i, j, k ) += -dt_ov_dz * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j, k-1 ) ) + dt_ov_dx * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i-1, j, k ) );
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i, j, k ) += -dt_ov_dx * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i-1, j, k ) ) + dt_ov_dy * ( ( *Ex3D )( i, j, k ) - ( *Ex3D )( i, j-1, k ) );
            }
        }
    }
    
}

