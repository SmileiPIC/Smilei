
#include "MF_Solver2D_Grassi.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MF_Solver2D_Grassi::MF_Solver2D_Grassi( Params &params )
    : Solver2D( params )
{

    dt_ov_dx  = dt/params.cell_length[0];
    dt_ov_dy  = dt/params.cell_length[0];
    dx        = params.cell_length[0];
    dy        = params.cell_length[1];
    if( dx!=dy ) {
        ERROR( "Grassi solver requires the same cell-length in x and y directions" );
    }
    
    //double alpha = 1.+(1.-(dt/dx)*(dt/dx))/8.;
    //double delta = (1.-alpha)/3.;
    double sigma = dt*dt/( dx*dx );
    double alpha = ( 5.-2.*sigma )/4.;
    double delta = ( 2.*sigma-1. )/12.;
    Ax    = alpha*dt/dx;
    Ay    = alpha*dt/dy;
    Dx    = delta*dt/dx;
    Dy    = delta*dt/dy;
    
    isEFilterApplied = params.Friedman_filter;
}

MF_Solver2D_Grassi::~MF_Solver2D_Grassi()
{
}

void MF_Solver2D_Grassi::operator()( ElectroMagn *fields )
{
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    // Static-cast of the fields
    Field2D *Ex2D;
    Field2D *Ey2D;
    Field2D *Ez2D;
    if( isEFilterApplied ) {
        Ex2D = static_cast<Field2D *>( fields->filter_->Ex_[0] );
        Ey2D = static_cast<Field2D *>( fields->filter_->Ey_[0] );
        Ez2D = static_cast<Field2D *>( fields->filter_->Ez_[0] );
    } else {
        Ex2D = static_cast<Field2D *>( fields->Ex_ );
        Ey2D = static_cast<Field2D *>( fields->Ey_ );
        Ez2D = static_cast<Field2D *>( fields->Ez_ );
    }
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    
    
    // Magnetic field Bx^(p,d) using Ez^(p,p)
    // --------------------------------------
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) { // j=0,1 & nx_d-2,nx_d-1 treated by exchange and/or BCs
        
            ( *Bx2D )( i, j ) += Ay * ( ( *Ez2D )( i, j-1 ) - ( *Ez2D )( i, j ) )
                                 +               Dy * ( ( *Ez2D )( i, j-2 ) - ( *Ez2D )( i, j+1 ) );
        }
    }
    
    
    // Magnetic field By^(d,p) using Ez^(p,p)
    // --------------------------------------
    for( unsigned int i=2 ; i<nx_d-2;  i++ ) { // i=0,1 & nx_d-2,nx_d-1 treated by exchange and/or BCs
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
        
            ( *By2D )( i, j ) += Ax * ( ( *Ez2D )( i, j )   - ( *Ez2D )( i-1, j ) )
                                 +               Dx * ( ( *Ez2D )( i+1, j ) - ( *Ez2D )( i-2, j ) );
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
    
    
    // Magnetic field Bz^(d,d) using Ex^(d,p) & Ey^(p,d)
    // -------------------------------------------------
    for( unsigned int i=2 ; i<nx_d-2;  i++ ) {      // i=0,1 & nx_d-2,nx_d-1 treated by exchange and/or BCs
        for( unsigned int j=2 ; j<ny_d-2 ; j++ ) {  // j=0,1 & nx_d-2,nx_d-1 treated by exchange and/or BCs
        
            ( *Bz2D )( i, j ) += Ax * ( ( *Ey2D )( i-1, j ) - ( *Ey2D )( i, j ) )
                                 +               Dx * ( ( *Ey2D )( i-2, j ) - ( *Ey2D )( i+1, j ) )
                                 +               Ay * ( ( *Ex2D )( i, j )   - ( *Ex2D )( i, j-1 ) )
                                 +               Dy * ( ( *Ex2D )( i, j+1 ) - ( *Ex2D )( i, j-2 ) );
                                 
        }
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
    
}

