
#include "MF_Solver2D_Cowan.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include <algorithm>

MF_Solver2D_Cowan::MF_Solver2D_Cowan( Params &params )
    : Solver2D( params )
{
    ERROR( "Under development, not yet working" );
    dx = params.cell_length[0];
    dy = params.cell_length[1];
    
    // Cowan parameters [see Nuter et al., Eur. Phys. J. D 68, 177 (2014) for details]
    double delta = std::min( dx, dy );
    double beta_xy =  0.125*delta*delta/( dy*dy );
    double beta_yx =  0.125*delta*delta/( dx*dx );
    double alpha_x =  1.-2.*beta_xy;
    double alpha_y =  1.-2.*beta_yx;
    
    Ax = dt_ov_dx * alpha_x;
    Ay = dt_ov_dy * alpha_y;
    Bx = dt_ov_dx * beta_xy;
    By = dt_ov_dy * beta_yx;
//    Ax = dt_ov_dx;
//    Ay = dt_ov_dy;
//    Bx = 0.;
//    By = 0.;

//    // Check for time filtering on Efields
//    istimeFilterApplied = false;
//    if (params.timeFilter_int!=0)
//        istimeFilterApplied = true;
}

MF_Solver2D_Cowan::~MF_Solver2D_Cowan()
{
}

void MF_Solver2D_Cowan::operator()( ElectroMagn *fields )
{
    // const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    // Static-cast of the fields
    Field2D *Ex2D;
    Field2D *Ey2D;
    Field2D *Ez2D;
//    if (istimeFilterApplied) {
//        Ex2D = static_cast<Field2D*>(fields->Ex_f);
//        Ey2D = static_cast<Field2D*>(fields->Ey_f);
//        Ez2D = static_cast<Field2D*>(fields->Ez_f);
//    } else {
    Ex2D = static_cast<Field2D *>( fields->Ex_ );
    Ey2D = static_cast<Field2D *>( fields->Ey_ );
    Ez2D = static_cast<Field2D *>( fields->Ez_ );
//    }
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    
    
    // Magnetic field Bx^(p,d)
    for( unsigned int i=1; i<nx_d-2;  i++ ) {
        for( unsigned int j=1; j<ny_d-1; j++ ) {
            ( *Bx2D )( i, j ) += Ay * ( ( *Ez2D )( i, j-1 ) - ( *Ez2D )( i, j ) )
                                 +               By * ( ( *Ez2D )( i+1, j-1 )-( *Ez2D )( i+1, j ) + ( *Ez2D )( i-1, j-1 )-( *Ez2D )( i-1, j ) );
        }//j
    }//i
    
    // Magnetic field By^(d,p) & Bz^(d,d)
    for( unsigned int i=1; i<nx_d-1; i++ ) {
    
        // By^(d,p)
        for( unsigned int j=1; j<ny_p-1; j++ ) {
            ( *By2D )( i, j ) += Ax * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i-1, j ) )
                                 +               Bx * ( ( *Ez2D )( i, j+1 )-( *Ez2D )( i-1, j+1 ) + ( *Ez2D )( i, j-1 )-( *Ez2D )( i-1, j-1 ) );
        }//j
        
        // Bz^(d,d)
        for( unsigned int j=1; j<ny_d-1; j++ ) {
            ( *Bz2D )( i, j ) += Ay * ( ( *Ex2D )( i, j ) - ( *Ex2D )( i, j-1 ) )
                                 +               By * ( ( *Ex2D )( i+1, j )-( *Ex2D )( i+1, j-1 ) + ( *Ex2D )( i-1, j )-( *Ex2D )( i-1, j-1 ) )
                                 -               Ax * ( ( *Ey2D )( i, j ) - ( *Ey2D )( i-1, j ) )
                                 -               Bx * ( ( *Ey2D )( i, j+1 )-( *Ey2D )( i-1, j+1 ) + ( *Ey2D )( i, j-1 )-( *Ey2D )( i-1, j-1 ) );
        }//j
    }//i
    
}//END solveMaxwellFaraday



