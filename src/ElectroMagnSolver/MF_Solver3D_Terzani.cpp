
#include "MF_Solver3D_Terzani.h"

#include "ElectroMagn.h"
#include "ElectroMagn3D.h"
#include "Field3D.h"

MF_Solver3D_Terzani::MF_Solver3D_Terzani( Params &params )
    : Solver3D( params )
{
    delta = ( dt_ov_dx * dt_ov_dx - 1.) / 12. ;
    // this coefficient and the associated solver are defined in 
    // D. Terzani and P. Londrillo, Computer Physics Communications 242 (2019)
    // https://doi.org/10.1016/j.cpc.2019.04.007
    // (Eq. 32 for the solver, coefficient delta defined after Eq. A5, choosing to modify only the MF solver)
}

MF_Solver3D_Terzani::~MF_Solver3D_Terzani()
{
}

void MF_Solver3D_Terzani::operator()( ElectroMagn *fields )
{
   const unsigned int nx_p = fields->dimPrim[0];
   const unsigned int nx_d = fields->dimDual[0];
   const unsigned int ny_p = fields->dimPrim[1];
   const unsigned int ny_d = fields->dimDual[1];
   const unsigned int nz_p = fields->dimPrim[2];
   const unsigned int nz_d = fields->dimDual[2];
  
  
    // Static-cast of the fields
    Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
    Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
    Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
    Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
    Field3D *By3D = static_cast<Field3D *>( fields->By_ );
    Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
    
    ElectroMagn3D *EM3D = static_cast<ElectroMagn3D *>( fields );
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *Bx3D )( i, j, k ) += -dt_ov_dy * ( ( *Ez3D )( i, j, k ) - ( *Ez3D )( i, j-1, k ) ) + dt_ov_dz * ( ( *Ey3D )( i, j, k ) - ( *Ey3D )( i, j, k-1 ) );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( i, j, k ) = ( *By3D )( i, j, k )
                                     -                 dt_ov_dz * ( ( *Ex3D )( i  , j, k ) - ( *Ex3D )( i  , j, k-1 ) ) 
                                     + (1.-3.*delta) * dt_ov_dx * ( ( *Ez3D )( i  , j, k ) - ( *Ez3D )( i-1, j, k   ) )
                                     + (      delta) * dt_ov_dx * ( ( *Ez3D )( i+1, j, k ) - ( *Ez3D )( i-2, j, k   ) );
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=2 ; i<nx_d-2 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( i, j, k ) = ( *Bz3D )( i, j, k )
                                     - (1.-3.*delta) * dt_ov_dx * ( ( *Ey3D )( i  , j, k ) - ( *Ey3D )( i-1, j  , k ) ) 
                                     - (      delta) * dt_ov_dx * ( ( *Ey3D )( i+1, j, k ) - ( *Ey3D )( i-2, j  , k ) ) 
                                     +                 dt_ov_dy * ( ( *Ex3D )( i  , j, k ) - ( *Ex3D )( i  , j-1, k ) );
            }
        }
    }
    
    // Left border: evolve as in Yee solver
    if( EM3D->isXmin ) {
    
        // Magnetic field By^(d,p,d)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( 1, j, k ) +=
                                      - dt_ov_dz * ( ( *Ex3D )( 1, j, k ) - ( *Ex3D )( 1, j, k-1 ) ) 
                                      + dt_ov_dx * ( ( *Ez3D )( 1, j, k ) - ( *Ez3D )( 0, j, k   ) );
            }
        }
        // Magnetic field Bz^(d,d,p)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( 1, j, k ) +=
                                      - dt_ov_dx * ( ( *Ey3D )( 1, j, k ) - ( *Ey3D )( 0, j  , k ) ) 
                                      + dt_ov_dy * ( ( *Ex3D )( 1, j, k ) - ( *Ex3D )( 1, j-1, k ) );
            }
        }

    } // end if isXmin
    
    // Right border: evolve as in Yee solver
    if( EM3D->isXmax ) {
    
        // Magnetic field By^(d,p,d)
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                ( *By3D )( nx_d-2, j, k ) +=
                                      - dt_ov_dz * ( ( *Ex3D )( nx_d-2, j, k ) - ( *Ex3D )( nx_d-2  , j  , k-1 ) ) 
                                      + dt_ov_dx * ( ( *Ez3D )( nx_d-2, j, k ) - ( *Ez3D )( nx_d-3, j  , k   ) );
            }
        }
        // Magnetic field Bz^(d,d,p)
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                ( *Bz3D )( nx_d-2, j, k ) +=
                                      - dt_ov_dx * ( ( *Ey3D )( nx_d-2, j, k ) - ( *Ey3D )( nx_d-3, j  , k ) ) 
                                      + dt_ov_dy * ( ( *Ex3D )( nx_d-2, j, k ) - ( *Ex3D )( nx_d-2  , j-1, k ) );
            }
        }

    } // end if isXmax
    
}

