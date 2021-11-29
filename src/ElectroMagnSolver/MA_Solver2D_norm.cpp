
#include "MA_Solver2D_norm.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MA_Solver2D_norm::MA_Solver2D_norm( Params &params )
    : Solver2D( params )
{
}

MA_Solver2D_norm::~MA_Solver2D_norm()
{
}

void MA_Solver2D_norm::operator()( ElectroMagn *fields )
{

    // Static-cast of the fields
    Field2D *Ex2D = static_cast<Field2D *>( fields->Ex_ );
    Field2D *Ey2D = static_cast<Field2D *>( fields->Ey_ );
    Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    Field2D *Jx2D = static_cast<Field2D *>( fields->Jx_ );
    Field2D *Jy2D = static_cast<Field2D *>( fields->Jy_ );
    Field2D *Jz2D = static_cast<Field2D *>( fields->Jz_ );

    // double sumJx = 0;
    // double sumJy = 0;
    // double sumJz = 0;

    // Electric field Ex^(d,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            ( *Ex2D )( i, j ) += -dt*( *Jx2D )( i, j ) + dt_ov_dy * ( ( *Bz2D )( i, j+1 ) - ( *Bz2D )( i, j ) );
            // sumJx += ( *Jx2D )( i, j );
        }
    }

    // Electric field Ey^(p,d)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            ( *Ey2D )( i, j ) += -dt*( *Jy2D )( i, j ) - dt_ov_dx * ( ( *Bz2D )( i+1, j ) - ( *Bz2D )( i, j ) );
            // sumJy += ( *Jy2D )( i, j );
        }
    }

    // Electric field Ez^(p,p)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            ( *Ez2D )( i, j ) += -dt*( *Jz2D )( i, j )
                                 +               dt_ov_dx * ( ( *By2D )( i+1, j ) - ( *By2D )( i, j ) )
                                 -               dt_ov_dy * ( ( *Bx2D )( i, j+1 ) - ( *Bx2D )( i, j ) );
            // sumJz += ( *Jz2D )( i, j );
        }
    }

    // std::cerr << std::scientific << std::setprecision( 10 )
    //           << "sum Jx: " << sumJx
    //           << " sum Jy: " << sumJy
    //           << " sum Jz: " << sumJz
    //           << std::endl;

}
