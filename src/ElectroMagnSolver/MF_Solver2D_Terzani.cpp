
#include "MF_Solver2D_Terzani.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MF_Solver2D_Terzani::MF_Solver2D_Terzani( Params &params )
    : Solver2D( params )
{
    isEFilterApplied = params.Friedman_filter;
    
    delta = ( dt_ov_dx * dt_ov_dx - 1.) / 12. ;
    // this coefficient and the associated solver are defined in 
    // D. Terzani and P. Londrillo, Computer Physics Communications 242 (2019)
    // https://doi.org/10.1016/j.cpc.2019.04.007
    // (Eq. 32 for the solver, coefficient delta defined after Eq. A5, choosing to modify only the MF solver)
}

MF_Solver2D_Terzani::~MF_Solver2D_Terzani()
{
}

void MF_Solver2D_Terzani::operator()( ElectroMagn *fields )
{
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    const unsigned int nz_p = fields->dimPrim[2];
    const unsigned int nz_d = fields->dimDual[2];
  
    // Static-cast of the fields
    Field2D *Ex2D;
    Field2D *Ey2D;
    if( !isEFilterApplied ) {
        Ex2D = static_cast<Field2D *>( fields->Ex_ );
        Ey2D = static_cast<Field2D *>( fields->Ey_ );
    } else {
        Ex2D = static_cast<Field2D *>( fields->filter_->Ex_[0] );
        Ey2D = static_cast<Field2D *>( fields->filter_->Ey_[0] );
    }
    Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
    
    // Magnetic field Bx^(p,d)
    for( unsigned int i=0 ; i<nx_d-1;  i++ ) {
        #pragma omp simd
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bx2D )( i, j ) -= dt_ov_dy * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i, j-1 ) );
        }
    }

    for( unsigned int i=2 ; i<nx_d-2;  i++ ) {
        // Magnetic field By^(d,p)
        #pragma omp simd
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            ( *By2D )( i, j ) = ( *By2D )( i, j )
                              + (1.-3.*delta) * dt_ov_dx * ( ( *Ez2D )( i  , j ) - ( *Ez2D )( i-1, j ) )
                              +       (delta) * dt_ov_dx * ( ( *Ez2D )( i+1, j ) - ( *Ez2D )( i-2, j ) );
        }
        
        // Magnetic field Bz^(d,d)
        #pragma omp simd
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            ( *Bz2D )( i, j ) = ( *Bz2D )( i, j )
                              +                 dt_ov_dy * ( ( *Ex2D )( i  , j ) - ( *Ex2D )( i  , j-1 ) )
                              - (1.-3.*delta) * dt_ov_dx * ( ( *Ey2D )( i  , j ) - ( *Ey2D )( i-1, j   ) )
                              -       (delta) * dt_ov_dx * ( ( *Ey2D )( i+1, j ) - ( *Ey2D )( i-2, j   ) );
        }
    }
    
    // Left border: evolve as in Yee solver
    
    unsigned int i = 1;
    // Magnetic field By^(d,p)
    #pragma omp simd
    for( unsigned int j=0 ; j<ny_p ; j++ ) {
        ( *By2D )( i, j ) += dt_ov_dx * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i-1, j ) );
    }
    // Magnetic field Bz^(d,d)
    #pragma omp simd
    for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
        ( *Bz2D )( i, j ) += dt_ov_dy * ( ( *Ex2D )( i, j ) - ( *Ex2D )( i  , j-1 ) )
                          -  dt_ov_dx * ( ( *Ey2D )( i, j ) - ( *Ey2D )( i-1, j   ) );
    }
    
    // Right border: evolve as in Yee solver
    i = nx_d-2;
    // Magnetic field By^(d,p)
    #pragma omp simd
    for( unsigned int j=0 ; j<ny_p ; j++ ) {
        ( *By2D )( i, j ) += dt_ov_dx * ( ( *Ez2D )( i, j ) - ( *Ez2D )( i-1, j ) ); 
    }
    // Magnetic field Bz^(d,d)
    #pragma omp simd
    for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
        ( *Bz2D )( i, j ) += dt_ov_dy * ( ( *Ex2D )( i, j ) - ( *Ex2D )( i  , j-1 ) )
                          -  dt_ov_dx * ( ( *Ey2D )( i, j ) - ( *Ey2D )( i-1, j   ) );
    }
    
}

