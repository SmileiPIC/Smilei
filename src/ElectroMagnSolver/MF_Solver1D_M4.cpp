
#include "MF_Solver1D_M4.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MF_Solver1D_M4::MF_Solver1D_M4( Params &params )
    : Solver1D( params )
{

    delta_x = ( dt_ov_dx * dt_ov_dx  - 1.0 ) / 12.0;
    alpha_x = 1.0 - 3.0 * delta_x;

    Ax    = alpha_x*dt_ov_dx;
    Dx    = delta_x*dt_ov_dx;
}

MF_Solver1D_M4::~MF_Solver1D_M4()
{
}

void MF_Solver1D_M4::operator()( ElectroMagn *fields )
{
    Field1D *Ey1D   = static_cast<Field1D *>( fields->Ey_ );
    Field1D *Ez1D   = static_cast<Field1D *>( fields->Ez_ );
    Field1D *By1D   = static_cast<Field1D *>( fields->By_ );
    Field1D *Bz1D   = static_cast<Field1D *>( fields->Bz_ );
    
    // ---------------------
    // Solve Maxwell-Faraday
    // ---------------------
    // NB: bx is given in 1d and defined when initializing the fields (here put to 0)
    // Transverse fields  by & bz are defined on the dual grid
    //for (unsigned int ix=1 ; ix<nx_p ; ix++) {
    for( unsigned int ix=2 ; ix<nx_d-2 ; ix++ ) {
        ( *By1D )( ix )= ( *By1D )( ix ) + Ax * ( ( *Ez1D )( ix   ) - ( *Ez1D )( ix-1 ) )
                                         + Dx * ( ( *Ez1D )( ix+1 ) - ( *Ez1D )( ix-2 ) ) ;
        ( *Bz1D )( ix )= ( *Bz1D )( ix ) - Ax * ( ( *Ey1D )( ix   ) - ( *Ey1D )( ix-1 ) )
                                         - Dx * ( ( *Ey1D )( ix+1 ) - ( *Ey1D )( ix-2 ) ) ;
    }
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    ( *By1D )( 1 )= ( *By1D )( 1 ) + dt_ov_dx * ( ( *Ez1D )( 1 ) - ( *Ez1D )( 0 ) );
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    ( *By1D )( nx_d-2 )= ( *By1D )( nx_d-2 ) + dt_ov_dx * ( ( *Ez1D )( nx_d-2 ) - ( *Ez1D )( nx_d-3 ) );
    // at Xmin+dx - treat using simple discretization of the curl (will be overwritten if not at the xmin-border)
    ( *Bz1D )( 1 )= ( *Bz1D )( 1 ) - dt_ov_dx * ( ( *Ey1D )( 1 ) - ( *Ey1D )( 0 ) );
    // at Xmax-dx - treat using simple discretization of the curl (will be overwritten if not at the xmax-border)
    ( *Bz1D )( nx_d-2 )= ( *Bz1D )( nx_d-2 ) - dt_ov_dx * ( ( *Ey1D )( nx_d-2 ) - ( *Ey1D )( nx_d-3 ) );
}

