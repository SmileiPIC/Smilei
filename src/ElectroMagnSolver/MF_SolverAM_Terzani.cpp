
#include "MF_SolverAM_Terzani.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverAM_Terzani::MF_SolverAM_Terzani( Params &params )
    : SolverAM( params )
{
    isEFilterApplied = params.Friedman_filter;
    
    delta = ( dt_ov_dl * dt_ov_dl - 1.) / 12. ;
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
        int  j_glob  = ( static_cast<ElectroMagnAM *>( fields ) )->j_glob_;
        bool isXmin  = ( static_cast<ElectroMagnAM *>( fields ) )->isXmin;
        bool isXmax  = ( static_cast<ElectroMagnAM *>( fields ) )->isXmax;
        bool isYmin  = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;

        // ---- Magnetic field Bl^(p,d)
        for( unsigned int i=0 ; i<nl_p;  i++ ) {
            #pragma omp simd
            for( unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bl )( i, j ) += - dt/( ( j_glob+j-0.5 )*dr ) * ( ( static_cast<double>( j+j_glob ))*( *Et )( i, j ) 
                                   - ( static_cast<double>(j+j_glob-1. ))*( *Et )( i, j-1 ) 
                                   + Icpx*(static_cast<double>(imode))   *( *Er )( i, j ) );
            }
        }

        // ---- Magnetic field Br^(d,p), smaller loop boundaries to avoid segfault
        for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * (1.-3.*delta) * ( ( *Et )( i  , j ) - ( *Et )( i-1, j ) )
                                 + dt_ov_dl * delta         * ( ( *Et )( i+1, j ) - ( *Et )( i-2, j ) )
                                +Icpx*dt*( static_cast<double>(imode) )/( ( static_cast<double>( j_glob+j ))*dr )*( *El )( i, j ) ;
            }
        }
        // Magnetic field Br^(d,p), left border of the patch: evolve as in Yee solver
        for( unsigned int i=1 ; i<2 ; i++ ) {
            #pragma omp simd 
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * ( ( *Et )( i, j ) - ( *Et )( i-1, j ) )
                                +Icpx*dt*( static_cast<double>(imode) )/( ( static_cast<double>( j_glob+j )) *dr )*( *El )( i, j ) ;
            }
        }
        // Magnetic field Br^(d,p), right border of the patch: evolve as in Yee solver
        for( unsigned int i=nl_d-2 ; i<nl_d-1 ; i++ ) {
            #pragma omp simd 
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) { //Specific condition on axis
                ( *Br )( i, j ) += dt_ov_dl * ( ( *Et )( i, j ) - ( *Et )( i-1, j ) )
                                 +Icpx*dt*( static_cast<double>(imode) )/( ( static_cast<double>( j_glob+j )) *dr )*( *El )( i, j ) ;
            }
        }
          
        // ---- Magnetic field Bt^(d,d), smaller loop boundaries to avoid segfault
        for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
            #pragma omp simd
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) += dt_ov_dr *                ( ( *El )( i, j   ) - ( *El )( i, j-1 ) )
                                  -dt_ov_dl * (1.-3*delta) * ( ( *Er )( i, j   ) - ( *Er )( i-1, j ) )
                                  -dt_ov_dl * delta        * ( ( *Er )( i+1, j ) - ( *Er )( i-2, j ) );
            }
        }
        
        // Magnetic field Bt^(d,d), left border of the patch : evolve as in Yee solver
        for( unsigned int i=1 ; i<2 ; i++ ) {
            #pragma omp simd 
            for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
                ( *Bt )( i, j ) += dt_ov_dr * ( ( *El )( i, j ) - ( *El )( i, j-1 ) )
                                  -dt_ov_dl * ( ( *Er )( i, j ) - ( *Er )( i-1, j ) );
            }
        }
        
        // Magnetic field Bt^(d,d), right border of the patch : evolve as in Yee solver
        for( unsigned int i=nl_d-2 ; i<nl_d-1 ; i++ ) {
          #pragma omp simd 
          for( unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++ ) {
              ( *Bt )( i, j ) +=  dt_ov_dr * ( ( *El )( i, j ) - ( *El )( i, j-1 ) )
                                 -dt_ov_dl * ( ( *Er )( i, j ) - ( *Er )( i-1, j ) );
          }
        }
        
        
        // On axis conditions: as in Yee solver, 
        // except for mode 1 where derivative along x is present
        if( isYmin ) {
            unsigned int j=2;
            if( imode==0 ) {
                // ---- Br
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    ( *Br )( i, j )= 0;
                    ( *Br )( i, 1 )=-( *Br )( i, 3 );
                }
                // ---- Bt
                for( unsigned int i= 0 ; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= -( *Bt )( i, j+1 );
                }
                // ---- Bl
                for( unsigned int i= 0 ; i<nl_p ; i++ ) {
                    ( *Bl )( i, j )= ( *Bl )( i, j+1 );
                }
            }

            else if( imode==1 ) {
                // ---- Bl
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    ( *Bl )( i, j )= -( *Bl )( i, j+1 ); // Zero Bl mode 1 on axis.
                }
                
                // ---- Br: use Terzani's derivative in smaller loop boundaries
                for( unsigned int i=2 ; i<nl_d-2 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                   +	 (1.-3.*delta)*	dt_ov_dl * ( ( *Et )( i  , j )-( *Et )( i-1, j ) )
                                   +   delta        * dt_ov_dl * ( ( *Et )( i+1, j )-( *Et )( i-2, j ) );
                    ( *Br )( i, 1 ) =( *Br )( i, 3 );
                }
                
                // ---- Br: as in Yee solver at the left border of the patch 
                for( unsigned int i=1 ; i<2 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                       +	dt_ov_dl*( ( *Et )( i, j   )-( *Et )( i-1, j ) );
                    ( *Br )( i, 1 ) =( *Br )( i, 3 );
                }
                
                // ---- Br: as in Yee solver at the right border of the patch 
                for( unsigned int i=nl_d-2 ; i<nl_d-1 ; i++ ) {
                    ( *Br )( i, j )+=  Icpx*dt_ov_dr*( *El )( i, j+1 )
                                       +	dt_ov_dl*( ( *Et )( i, j   )-( *Et )( i-1, j ) );
                    ( *Br )( i, 1 ) =( *Br )( i, 3 );
                }
                
                // ---- Bt
                for( unsigned int i=0; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= ( *Bt )( i, j+1 ); // Non zero Bt mode 1 on axis.
                }

            } else { // modes > 1
                // ---- Bl
                for( unsigned int  i=0 ; i<nl_p; i++ ) {
                    ( *Bl )( i, j )= -( *Bl )( i, j+1 );
                }
                // ---- Br
                for( unsigned int i=0 ; i<nl_d; i++ ) {
                    ( *Br )( i, j )= 0;
                    ( *Br )( i, 1 )=-( *Br )( i, 3 );
                }
                // ---- Bt
                for( unsigned int  i=0 ; i<nl_d ; i++ ) {
                    ( *Bt )( i, j )= - ( *Bt )( i, j+1 );
                }
            }
        }
    }
}
