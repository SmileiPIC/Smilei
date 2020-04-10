
#include "MA_SolverAM_norm.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"
#include "Patch.h"
MA_SolverAM_norm::MA_SolverAM_norm( Params &params )
    : SolverAM( params )
{
}

MA_SolverAM_norm::~MA_SolverAM_norm()
{
}

void MA_SolverAM_norm::operator()( ElectroMagn *fields )
{
    for( unsigned int imode=0 ; imode<Nmode ; imode++ ) {
    
        // Static-cast of the fields_SolverAM_norm.cpp
        cField2D *El = ( static_cast<ElectroMagnAM *>( fields ) )->El_[imode];
        cField2D *Er = ( static_cast<ElectroMagnAM *>( fields ) )->Er_[imode];
        cField2D *Et = ( static_cast<ElectroMagnAM *>( fields ) )->Et_[imode];
        cField2D *Bl = ( static_cast<ElectroMagnAM *>( fields ) )->Bl_[imode];
        cField2D *Br = ( static_cast<ElectroMagnAM *>( fields ) )->Br_[imode];
        cField2D *Bt = ( static_cast<ElectroMagnAM *>( fields ) )->Bt_[imode];
        cField2D *Jl = ( static_cast<ElectroMagnAM *>( fields ) )->Jl_[imode];
        cField2D *Jr = ( static_cast<ElectroMagnAM *>( fields ) )->Jr_[imode];
        cField2D *Jt = ( static_cast<ElectroMagnAM *>( fields ) )->Jt_[imode];
        int j_glob    = ( static_cast<ElectroMagnAM *>( fields ) )->j_glob_;
        bool isYmin = ( static_cast<ElectroMagnAM *>( fields ) )->isYmin;
        //double *invR = ( static_cast<ElectroMagnAM *>( fields ) )->invR;
        //double *invRd = ( static_cast<ElectroMagnAM *>( fields ) )->invRd;
        
        // Electric field Elr^(d,p)
        for( unsigned int i=0 ; i<nl_d ; i++ ) {
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
                ( *El )( i, j ) += -dt*( *Jl )( i, j )
                                   +                 dt/( ( j_glob+j )*dr )*( ( j+j_glob+0.5 )*( *Bt )( i, j+1 ) - ( j+j_glob-0.5 )*( *Bt )( i, j ) )
                                   +                 Icpx*dt*( double )imode/( ( j_glob+j )*dr )*( *Br )( i, j );
            }
        }
        for( unsigned int i=0 ; i<nl_p ; i++ ) {
            for( unsigned int j=isYmin*3 ; j<nr_d ; j++ ) {
                ( *Er )( i, j ) += -dt*( *Jr )( i, j )
                                   -                  dt_ov_dl * ( ( *Bt )( i+1, j ) - ( *Bt )( i, j ) )
                                   -                  Icpx*dt*( double )imode/( ( j_glob+j-0.5 )*dr )* ( *Bl )( i, j );
                                   
            }
        }
        for( unsigned int i=0 ;  i<nl_p ; i++ ) {
            for( unsigned int j=isYmin*3 ; j<nr_p ; j++ ) {
                ( *Et )( i, j ) += -dt*( *Jt )( i, j )
                                   +                  dt_ov_dl * ( ( *Br )( i+1, j ) - ( *Br )( i, j ) )
                                   -                  dt_ov_dr * ( ( *Bl )( i, j+1 ) - ( *Bl )( i, j ) );
            }
        }
        if( isYmin ) { 
            // Conditions on axis
            unsigned int j=2;
            if( imode==0 ) {
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    ( *Et )( i, j )=0;
                    ( *Et )( i, j-1 )=-( *Et )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    //( *Er )( i, j+1 )= ( *Er )( i, j+2 ) / 9.;
                    ( *Er )( i, j )= -( *Er )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_d ; i++ ) {
                    ( *El )( i, j )+= 4.*dt_ov_dr*( *Bt )( i, j+1 )-dt*( *Jl )( i, j );
                    ( *El )( i, j-1 )=( *El )( i, j+1 );
                }
            } else if( imode==1 ) {
                for( unsigned int i=0 ; i<nl_d  ; i++ ) {
                    ( *El )( i, j )= 0;
                    ( *El )( i, j-1 )=-( *El )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p  ; i++ ) {
                    //( *Et )( i, j )= -1./3.*( 4.*Icpx*( *Er )( i, j+1 )+( *Et )( i, j+1 ) );
                    ( *Et )( i, j )= -Icpx/8.*( 9.*( *Er )( i, j+1 )-( *Er )( i, j+2 ) );
                    ( *Et )( i, j-1 )=( *Et )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p ; i++ ) {
                    ( *Er )( i, j )=2.*Icpx*( *Et )( i, j )-( *Er )( i, j+1 );
                }
            } else { // mode > 1
                for( unsigned int  i=0 ; i<nl_d; i++ ) {
                    ( *El )( i, j )= 0;
                    ( *El )( i, j-1 )=-( *El )( i, j+1 );
                }
                for( unsigned int  i=0 ; i<nl_p; i++ ) {
                    ( *Er )( i, j+1 )= ( *Er )( i, j+2 ) / 9.;
                    ( *Er )( i, j )= -( *Er )( i, j+1 );
                }
                for( unsigned int i=0 ; i<nl_p; i++ ) {
                    ( *Et )( i, j )= 0;
                    ( *Et )( i, j-1 )=-( *Et )( i, j+1 );
                }
            }
        }
    }
}

