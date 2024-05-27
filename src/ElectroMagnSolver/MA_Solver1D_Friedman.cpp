
#include "MA_Solver1D_Friedman.h"

#include "ElectroMagn.h"
#include "Field1D.h"

MA_Solver1D_Friedman::MA_Solver1D_Friedman( Params &params )
    : Solver1D( params )
{
    ftheta = params.Friedman_theta;
    alpha  = 1.-0.5*ftheta+0.5*ftheta*ftheta;
    beta   = ftheta*( 1.-0.5*ftheta );
    delta  = 0.5*ftheta*( 1.-ftheta )*( 1.-ftheta );

    //std::cout << " FRIEDMAN FILTER PARAMETER : " << ftheta << std::endl;
}

MA_Solver1D_Friedman::~MA_Solver1D_Friedman()
{
}

void MA_Solver1D_Friedman::operator()( ElectroMagn *fields )
{

    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    // Static-cast of the fields
    Field1D *Ex1D   = static_cast<Field1D *>( fields->Ex_ );
    Field1D *Ey1D   = static_cast<Field1D *>( fields->Ey_ );
    Field1D *Ez1D   = static_cast<Field1D *>( fields->Ez_ );
    // Field1D *Bx1D   = static_cast<Field1D *>( fields->Bx_ );
    Field1D *By1D   = static_cast<Field1D *>( fields->By_ );
    Field1D *Bz1D   = static_cast<Field1D *>( fields->Bz_ );
    Field1D *Jx1D   = static_cast<Field1D *>( fields->Jx_ );
    Field1D *Jy1D   = static_cast<Field1D *>( fields->Jy_ );
    Field1D *Jz1D   = static_cast<Field1D *>( fields->Jz_ );

    Field1D *Ex_f = static_cast<Field1D *>( fields->filter_->Ex_[0] );
    Field1D *Ey_f = static_cast<Field1D *>( fields->filter_->Ey_[0] );
    Field1D *Ez_f = static_cast<Field1D *>( fields->filter_->Ez_[0] );
    Field1D *Ex_m1 = static_cast<Field1D *>( fields->filter_->Ex_[1] );
    Field1D *Ey_m1 = static_cast<Field1D *>( fields->filter_->Ey_[1] );
    Field1D *Ez_m1 = static_cast<Field1D *>( fields->filter_->Ez_[1] );
    Field1D *Ex_m2 = static_cast<Field1D *>( fields->filter_->Ex_[2] );
    Field1D *Ey_m2 = static_cast<Field1D *>( fields->filter_->Ey_[2] );
    Field1D *Ez_m2 = static_cast<Field1D *>( fields->filter_->Ez_[2] );

    double adv = 0.;

    // Electric field Ex^(d,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {

            adv               = -dt*( *Jx1D )( i ) ;
            // advance electric field
            ( *Ex1D )( i )   += adv;
            // compute the time-filtered field
            ( *Ex_f )( i )    = alpha*( *Ex1D )( i ) + beta*adv + delta*( ( *Ex_m1 )( i )+ftheta*( *Ex_m2 )( i ) );
            // update Ex_m2 and Ex_m1
            ( *Ex_m2 )( i )   = ( *Ex_m1 )( i ) - ftheta*( *Ex_m2 )( i );
            ( *Ex_m1 )( i )   = ( *Ex1D )( i )  - adv;

    }


    // Electric field Ey^(p,d)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {

            adv               = -dt*( *Jy1D )( i ) - dt_ov_dx * ( ( *Bz1D )( i+1 ) - ( *Bz1D )( i ) );
            // advance electric field
            ( *Ey1D )( i )   += adv;
            // compute the time-filtered field
            ( *Ey_f )( i )    = alpha*( *Ey1D )( i ) + beta*adv + delta*( ( *Ey_m1 )( i )+ftheta*( *Ey_m2 )( i ) );
            // update Ex_m2 and Ex_m1
            ( *Ey_m2 )( i )   = ( *Ey_m1 )( i ) - ftheta*( *Ey_m2 )( i ) ;
            ( *Ey_m1 )( i )   = ( *Ey1D )( i )  - adv;

    }


    // Electric field Ez^(p,p)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {

            adv               = - dt*( *Jz1D )( i )
                              + dt_ov_dx * ( ( *By1D )( i+1 ) - ( *By1D )( i ) );
            // advance electric field
            ( *Ez1D )( i )   += adv;
            // compute the time-filtered field
            ( *Ez_f )( i )    = alpha*( *Ez1D )( i ) + beta*adv + delta*( ( *Ez_m1 )( i )+ftheta*( *Ez_m2 )( i ) );
            // update Ex_m2 and Ex_m1
            ( *Ez_m2 )( i )   = ( *Ez_m1 )( i ) - ftheta * ( *Ez_m2 )( i ) ;
            ( *Ez_m1 )( i )   = ( *Ez1D )( i ) - adv ;

    }

}
