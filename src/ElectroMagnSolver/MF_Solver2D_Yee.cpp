
#include "MF_Solver2D_Yee.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MF_Solver2D_Yee::MF_Solver2D_Yee( Params &params )
    : Solver2D( params )
{
    isEFilterApplied = false;
    if( params.Friedman_filter ) {
        isEFilterApplied = true;
    }
}

MF_Solver2D_Yee::~MF_Solver2D_Yee()
{
}

void MF_Solver2D_Yee::operator()( ElectroMagn *fields )
{
    // Static-cast of the fields
    Field2D *Ex2D;
    Field2D *Ey2D;
    if( !isEFilterApplied ) {
        Ex2D = static_cast<Field2D *>( fields->Ex_ );
        Ey2D = static_cast<Field2D *>( fields->Ey_ );
    } else {
        Ex2D = static_cast<Field2D *>( fields->Exfilter[0] );
        Ey2D = static_cast<Field2D *>( fields->Eyfilter[0] );
    }
    Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    
    // Magnetic field Bx^(p,d)
    //cout << "nx_p,nx_d-1" << nx_p << " " << nx_d-1 ;
    //{
    {
        #pragma omp simd
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bx2D )( 0, j ) -= dt_ov_dy * ( ( *Ez2D )( 0, j ) - ( *Ez2D )( 0, j-1 ) );
        }
    }
    //    for (unsigned int i=0 ; i<nx_p;  i++) {
    for( unsigned int i=1 ; i<nx_d-1;  i++ ) {
        #pragma omp simd
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bx2D )( i, j ) -= dt_ov_dy * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i, j-1 ) );
        }
        //    }
        
        // Magnetic field By^(d,p)
        //    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            ( *By2D )( i, j ) += dt_ov_dx * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i-1, j ) );
        }
        //}
        
        // Magnetic field Bz^(d,d)
        //for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bz2D )( i, j ) += dt_ov_dy * ( ( *Ex2D )( i, j ) - ( *Ex2D )( i, j-1 ) )
                                 -               dt_ov_dx * ( ( *Ey2D )( i, j ) - ( *Ey2D )( i-1, j ) );
        }
    }
    //}// end parallel
}

