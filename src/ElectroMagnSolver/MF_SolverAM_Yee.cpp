
#include "MF_SolverAM_Yee.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverAM_Yee::MF_SolverAM_Yee( Params &params )
    : SolverAM( params )
{
    isEFilterApplied = false;
    if( params.Friedman_filter ) {
        ERROR( "Filters are not available yet" );
    };
}

MF_SolverAM_Yee::~MF_SolverAM_Yee()
{
}

void MF_SolverAM_Yee::operator()( ElectroMagn *fields )
{

    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
    
        // Static-cast of the fields
        
        cField2D *El = ( static_cast<ElectroMagnAM *>( fields ) )->El_[imode];
        cField2D *Er = ( static_cast<ElectroMagnAM *>( fields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( fields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( fields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( fields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( fields ) )->Bt_[imode];
        int  j_glob = ( static_cast<ElectroMagnAM *>( fields ) )->j_glob_;
        bool isYmin = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;
        //double *invR = ( static_cast<ElectroMagnAM *>( fields ) )->invR;
        //double *invRd = ( static_cast<ElectroMagnAM *>( fields ) )->invRd;
        
        // Magnetic field Bl^(p,d)
        for( unsigned int i=0 ; i<nl_p;  i++ ) {
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += - dt/( ( j_glob+j-0.5 )*dr ) * ( ( double )( j+j_glob )*( *Et )( i, j ) - ( double )( j+j_glob-1. )*( *Et )( i, j-1 ) + Icpx*( double )imode*( *Er )( i, j ) );
            }
        }
        
        // Magnetic field Br^(d,p)
        for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * ( ( *Et )( i, j ) - ( *Et )( i-1, j ) )
                                   +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i, j ) ;
            }
        }
        // Magnetic field Bt^(d,d)
        for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) += dt_ov_dr * ( ( *El )( i, j ) - ( *El )( i, j-1 ) )
                                   -dt_ov_dl * ( ( *Er )( i, j ) - ( *Er )( i-1, j ) );
            }
        }
        
        // On axis conditions
        if( isYmin ) {
            unsigned int j=2;
            if( imode==0 ) {
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    ( *Br )( i, j )=0;
                    ( *Br )( i, 1 )=-( *Br )( i, 3 );
                }
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    //( *Bt )( i, j+1 )= ( *Bt )( i, j+2 )/9.;
                    ( *Bt )( i, j )= -( *Bt )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p ; i++ ) {
                    ( *Bl )( i, j )= ( *Bl )( i, j+1 );
                }
            }
            
            else if( imode==1 ) {
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    ( *Bl )( i, j )= -( *Bl )( i, j+1 );
                }
                
                for( unsigned int i=1 ; i<nl_d-1 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                       +			dt_ov_dl*( ( *Et )( i, j )-( *Et )( i-1, j ) );
                    ( *Br )( i, 1 )=( *Br )( i, 3 );
                }
                for( unsigned int i=0; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= -2.*Icpx*( *Br )( i, j )-( *Bt )( i, j+1 );
                }
                
            } else { // modes > 1
                for( unsigned int  i=0 ; i<nl_p; i++ ) {
                    ( *Bl )( i, j )= -( *Bl )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_d; i++ ) {
                    ( *Br )( i, j )= 0;
                    ( *Br )( i, 1 )=-( *Br )( i, 3 );
                }
                for( unsigned int  i=0 ; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= - ( *Bt )( i, j+1 );
                }
            }
        }
    }
}

