
#include "MF_SolverAM_Lehe.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverAM_Lehe::MF_SolverAM_Lehe( Params &params )
    : SolverAM( params )
{

    beta_rl = 0.25; 
    beta_tl = 0.25;
    delta_l = 0.25 * ( 1.- ( sin( M_PI * dt_ov_dl * 0.5 )/dt_ov_dl )*( sin( M_PI*dt_ov_dl * 0.5 )/dt_ov_dl ) );

    alpha_r =  1. - 2.*beta_rl; 
    alpha_t =  1. - 2.*beta_tl; 
    alpha_l =  1. - 3.*delta_l;
    
}

MF_SolverAM_Lehe::~MF_SolverAM_Lehe()
{
}

void MF_SolverAM_Lehe::operator()( ElectroMagn *fields )
{

    const unsigned int nl_p = fields->dimPrim[0];
    const unsigned int nl_d = fields->dimDual[0];
    const unsigned int nr_p = fields->dimPrim[1];
    const unsigned int nr_d = fields->dimDual[1];
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
        bool isXmin = ( static_cast<ElectroMagnAM *>( fields ) )->isXmin;
        bool isXmax = ( static_cast<ElectroMagnAM *>( fields ) )->isXmax;
        //double *invR = ( static_cast<ElectroMagnAM *>( fields ) )->invR;
        //double *invRd = ( static_cast<ElectroMagnAM *>( fields ) )->invRd;

        // Magnetic field Bl^(p,d)
        for( unsigned int i=1 ; i<nl_p-1;  i++ ) {
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += -dt/( ( j_glob+j-0.5 )*dr ) * alpha_r * ( ( double )( j+j_glob )*( *Et )( i  , j ) - ( double )( j+j_glob-1. )*( *Et )( i  , j-1 ) )
                                   -dt/( ( j_glob+j-0.5 )*dr ) * beta_rl * ( ( double )( j+j_glob )*( *Et )( i+1, j ) - ( double )( j+j_glob-1. )*( *Et )( i+1, j-1 ) )
                                   -dt/( ( j_glob+j-0.5 )*dr ) * beta_rl * ( ( double )( j+j_glob )*( *Et )( i-1, j ) - ( double )( j+j_glob-1. )*( *Et )( i-1, j-1 ) )         
                                   -dt/( ( j_glob+j-0.5 )*dr ) * alpha_t * Icpx * ( double )imode  *( *Er )( i  , j ) 
                                   -dt/( ( j_glob+j-0.5 )*dr ) * beta_tl * Icpx * ( double )imode  *( *Er )( i+1, j )                  
                                   -dt/( ( j_glob+j-0.5 )*dr ) * beta_tl * Icpx * ( double )imode  *( *Er )( i-1, j ) ;
            }
        }

        // Magnetic field Br^(d,p)
        for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) +=  dt_ov_dl * alpha_l * ( ( *Et )( i  , j ) - ( *Et )( i-1, j ) )
                                   +dt_ov_dl * delta_l * ( ( *Et )( i+1, j ) - ( *Et )( i-2, j ) )
                                   +Icpx*dt*alpha_t*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i  , j )
                                   +Icpx*dt*beta_tl*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i+1, j ) 
                                   +Icpx*dt*beta_tl*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i-1, j ) ; 
            }
        }
        // Magnetic field Bt^(d,d)
        for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) +=  dt_ov_dr * alpha_r * ( ( *El )( i  , j ) - ( *El )( i  , j-1 ) )
                                   +dt_ov_dr * beta_rl * ( ( *El )( i+1, j ) - ( *El )( i+1, j-1 ) )
                                   +dt_ov_dr * beta_rl * ( ( *El )( i-1, j ) - ( *El )( i-1, j-1 ) )
                                   -dt_ov_dl * alpha_l * ( ( *Er )( i  , j ) - ( *Er )( i-1, j   ) )
                                   -dt_ov_dl * delta_l * ( ( *Er )( i+1, j ) - ( *Er )( i-2, j   ) );
            }
        }

        if ( isXmin ) { // at the left border, use a the classic finite differences of Yee solver
            unsigned int i=0;
            // Bl
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += -dt/( ( j_glob+j-0.5 )*dr ) * ( ( double )( j+j_glob )*( *Et )( i  , j ) - ( double )( j+j_glob-1. )*( *Et )( i  , j-1 ) )
                                   -dt/( ( j_glob+j-0.5 )*dr ) * Icpx * ( double )imode  *( *Er )( i  , j );
            }
            // Br
            i = 1;
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) +=  dt_ov_dl * ( ( *Et )( i  , j ) - ( *Et )( i-1, j ) )
                                   +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i  , j ); 
            } 
            // Bt
            i = 1;
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) +=  dt_ov_dr * ( ( *El )( i  , j ) - ( *El )( i  , j-1 ) )
                                   -dt_ov_dl * ( ( *Er )( i  , j ) - ( *Er )( i-1, j   ) );
            }
          
        }
        
        if ( isXmax ) { // at the right border, use a the classic finite differences of Yee solver
            
            // Bl
            unsigned int i=nl_p-1;
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += -dt/( ( j_glob+j-0.5 )*dr ) * ( ( double )( j+j_glob )*( *Et )( i  , j ) - ( double )( j+j_glob-1. )*( *Et )( i  , j-1 ) )
                                   -dt/( ( j_glob+j-0.5 )*dr ) * Icpx * ( double )imode  *( *Er )( i  , j );
            }
            // Br
            i=nl_d-2;
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { 
                ( *Br )( i, j ) +=  dt_ov_dl * ( ( *Et )( i  , j ) - ( *Et )( i-1, j ) )
                                   +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i  , j ); 
            } 
            // Bt
            i=nl_d-2;
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) +=  dt_ov_dr * ( ( *El )( i  , j ) - ( *El )( i  , j-1 ) )
                                   -dt_ov_dl * ( ( *Er )( i  , j ) - ( *Er )( i-1, j   ) );
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

                for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                       + dt_ov_dl * alpha_l* ( ( *Et )( i, j )-( *Et )( i-1, j ) )
                                       + dt_ov_dl * delta_l* ( ( *Et )( i+1, j ) - ( *Et )( i-2, j ) );
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
