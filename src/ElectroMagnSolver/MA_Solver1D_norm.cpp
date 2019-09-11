
#include "MA_Solver1D_norm.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MA_Solver1D_norm::MA_Solver1D_norm( Params &params )
    : Solver1D( params )
{
}

MA_Solver1D_norm::~MA_Solver1D_norm()
{
}

void MA_Solver1D_norm::operator()( ElectroMagn *fields )
{
    Field1D *Ex1D = static_cast<Field1D *>( fields->Ex_ );
    Field1D *Ey1D = static_cast<Field1D *>( fields->Ey_ );
    Field1D *Ez1D = static_cast<Field1D *>( fields->Ez_ );
    Field1D *By1D = static_cast<Field1D *>( fields->By_ );
    Field1D *Bz1D = static_cast<Field1D *>( fields->Bz_ );
    Field1D *Jx1D = static_cast<Field1D *>( fields->Jx_ );
    Field1D *Jy1D = static_cast<Field1D *>( fields->Jy_ );
    Field1D *Jz1D = static_cast<Field1D *>( fields->Jz_ );
    
    // --------------------
    // Solve Maxwell-Ampere
    // --------------------
    // Calculate the electrostatic field ex on the dual grid
    for( unsigned int ix=0 ; ix<nx_d ; ix++ ) {
        ( *Ex1D )( ix )= ( *Ex1D )( ix ) - dt * ( *Jx1D )( ix ) ;
    }
    // Transverse fields ey, ez  are defined on the primal grid
    for( unsigned int ix=0 ; ix<nx_p ; ix++ ) {
        ( *Ey1D )( ix )= ( *Ey1D )( ix ) - dt_ov_dx * ( ( *Bz1D )( ix+1 ) - ( *Bz1D )( ix ) ) - dt * ( *Jy1D )( ix ) ;
        ( *Ez1D )( ix )= ( *Ez1D )( ix ) + dt_ov_dx * ( ( *By1D )( ix+1 ) - ( *By1D )( ix ) ) - dt * ( *Jz1D )( ix ) ;
    }
}

