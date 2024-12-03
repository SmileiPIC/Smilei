
#include "MF_SolverAM_Terzani.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverAM_Terzani::MF_SolverAM_Terzani( Params &params )
    : SolverAM( params )
{
    isEFilterApplied = params.Friedman_filter;
    
    delta = ( pow(dt_ov_dl,2)-1.) / 12. ;
    // this coefficient and the associated solver are defined in 
    // D. Terzani and P. Londrillo, Computer Physics Communications 242 (2019)
    // https://doi.org/10.1016/j.cpc.2019.04.007
    // (Eq. 32 for the solver, coefficient delta defined after Eq. A5, choosing to modify only the MF solver)
}

MF_SolverAM_Terzani::~MF_SolverAM_Terzani()
{
}

void MF_SolverAM_Terzani::operator()( ElectroMagn *fields )
{

    const unsigned int nl_p = fields->dimPrim[0];
    const unsigned int nl_d = fields->dimDual[0];
    const unsigned int nr_p = fields->dimPrim[1];
    const unsigned int nr_d = fields->dimDual[1];
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {

        // Static-cast of the fields

        cField2D *El;
        cField2D *Er;
        cField2D *Et;
        if (isEFilterApplied) {
            El = static_cast<cField2D*>(fields->filter_->El_[imode][0]);
            Er = static_cast<cField2D*>(fields->filter_->Er_[imode][0]);
            Et = static_cast<cField2D*>(fields->filter_->Et_[imode][0]);
        } else {
            El = ( static_cast<ElectroMagnAM *>( fields ) )->El_[imode];
            Er = ( static_cast<ElectroMagnAM *>( fields ) )->Er_[imode];
            Et = ( static_cast<ElectroMagnAM *>( fields ) )->Et_[imode];
        }
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( fields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( fields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( fields ) )->Bt_[imode];
        int  j_glob = ( static_cast<ElectroMagnAM *>( fields ) )->j_glob_;
        bool isXmin = ( static_cast<ElectroMagnAM *>( fields ) )->isXmin;
        bool isXmax = ( static_cast<ElectroMagnAM *>( fields ) )->isXmax;
        bool isYmin = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;

        // Magnetic field Bl^(p,d)
        for( unsigned int i=0 ; i<nl_p;  i++ ) {
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += - dt/( ( j_glob+j-0.5 )*dr ) * ( ( double )( j+j_glob )*( *Et )( i, j ) - ( double )( j+j_glob-1. )*( *Et )( i, j-1 ) + Icpx*( double )imode*( *Er )( i, j ) );
            }
        }

        // Magnetic field Br^(d,p)
        for( unsigned int i=2+1*isXmin ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * (1.-3.*delta) * ( ( *Et )( i  , j ) - ( *Et )( i-1, j ) )
                                 + dt_ov_dl * delta        * ( ( *Et )( i+1, j ) - ( *Et )( i-2, j ) )
                                +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i, j ) ;
            }
        }
        
        if (isXmin){
            // Magnetic field Br^(d,p), left border: evolve as in Yee solver
          for( unsigned int i=1 ; i<3 ; i++ ) {
            #pragma omp simd 
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * ( ( *Et )( i, j ) - ( *Et )( i-1, j ) )
                                 +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i, j ) ;
            }
          }
        }
        
        // if (isXmax){
        //     // Magnetic field Br^(d,p), right border: evolve as in Yee solver
        //     unsigned int i=nl_d-2;
        //     #pragma omp simd 
        //     for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
        //         ( *Br )( i, j ) += dt_ov_dl * ( ( *Et )( i, j ) - ( *Et )( i-1, j ) )
        //                          +Icpx*dt*( double )imode/( ( double )( j_glob+j )*dr )*( *El )( i, j ) ;
        //     }
        // }
        
        
        // Magnetic field Bt^(d,d)
        for( unsigned int i=2+1*isXmin ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) += dt_ov_dr *                ( ( *El )( i, j   ) - ( *El )( i, j-1 ) )
                                  -dt_ov_dl * (1.-3*delta) * ( ( *Er )( i, j   ) - ( *Er )( i-1, j ) )
                                  -dt_ov_dl * delta        * ( ( *Er )( i+1, j ) - ( *Er )( i-2, j ) );
            }
        }
        
        if (isXmin){
            // Magnetic field Bt^(d,d), left border: evolve as in Yee solver
          for( unsigned int i=1 ; i<3 ; i++ ) {
            #pragma omp simd 
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) += dt_ov_dr * ( ( *El )( i, j ) - ( *El )( i, j-1 ) )
                                   -dt_ov_dl* ( ( *Er )( i, j ) - ( *Er )( i-1, j ) );
            }
          }
        }
        // if (isXmax){
        //     // Magnetic field Bt^(d,d), right border: evolve as in Yee solver
        //     unsigned int i=nl_d-2;
        //     #pragma omp simd 
        //     for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
        //         ( *Bt )( i, j ) += dt_ov_dr * ( ( *El )( i, j ) - ( *El )( i, j-1 ) )
        //                            -dt_ov_dl* ( ( *Er )( i, j ) - ( *Er )( i-1, j ) );
        //     }
        // }

        // On axis conditions
        if( isYmin ) {
            unsigned int j=2;
            if( imode==0 ) {
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    ( *Br )( i, j )=0;
                    ( *Br )( i, 1 )=-( *Br )( i, 3 );
                }
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= -( *Bt )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p ; i++ ) {
                    ( *Bl )( i, j )= ( *Bl )( i, j+1 );
                }
            }

            else if( imode==1 ) {
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    ( *Bl )( i, j )= -( *Bl )( i, j+1 ); // Zero Bl mode 1 on axis.
                }

                for( unsigned int i=2+1*isXmin ; i<nl_d-2 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                    +	 (1.-3.*delta) *	dt_ov_dl*( ( *Et )( i, j   )-( *Et )( i-1, j ) )
                                    +   delta       * dt_ov_dl*( ( *Et )( i+1, j )-( *Et )( i-2, j ) );
                    ( *Br )( i, 1 )=( *Br )( i, 3 );
                }
                
                if (isXmin){
                    // left border: as in Yee solver
                  for( unsigned int i=1 ; i<3 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                       +	dt_ov_dl*( ( *Et )( i, j   )-( *Et )( i-1, j ) );
                    ( *Br )( i, 1 )=( *Br )( i, 3 );
                  }
                }
                // if (isXmax){        
                //     // right border
                //     ( *Br )( nl_d-2, j )+=  Icpx*dt_ov_dr*( *El )( nl_d-2, j+1 )
                //                         +	dt_ov_dl*( ( *Et )( nl_d-2, j   )-( *Et )( nl_d-3, j ) );
                //     ( *Br )( nl_d-2, 1 )=( *Br )( nl_d-2, 3 );
                // } 
                
                for( unsigned int i=0; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= ( *Bt )( i, j+1 ); // Non zero Bt mode 1 on axis.
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
