
#include "MA_Solver2D_Friedman.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MA_Solver2D_Friedman::MA_Solver2D_Friedman( Params &params )
    : Solver2D( params )
{
    ftheta = params.Friedman_theta;
    alpha  = 1.-0.5*ftheta+0.5*ftheta*ftheta;
    beta   = ftheta*( 1.-0.5*ftheta );
    delta  = 0.5*ftheta*( 1.-ftheta )*( 1.-ftheta );
    
    //std::cout << " FRIEDMAN FILTER PARAMETER : " << ftheta << std::endl;
}

MA_Solver2D_Friedman::~MA_Solver2D_Friedman()
{
}

void MA_Solver2D_Friedman::operator()( ElectroMagn *fields )
{

    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    // Static-cast of the fields
    Field2D *Ex2D   = static_cast<Field2D *>( fields->Ex_ );
    Field2D *Ey2D   = static_cast<Field2D *>( fields->Ey_ );
    Field2D *Ez2D   = static_cast<Field2D *>( fields->Ez_ );
    Field2D *Bx2D   = static_cast<Field2D *>( fields->Bx_ );
    Field2D *By2D   = static_cast<Field2D *>( fields->By_ );
    Field2D *Bz2D   = static_cast<Field2D *>( fields->Bz_ );
    Field2D *Jx2D   = static_cast<Field2D *>( fields->Jx_ );
    Field2D *Jy2D   = static_cast<Field2D *>( fields->Jy_ );
    Field2D *Jz2D   = static_cast<Field2D *>( fields->Jz_ );
    
    Field2D *Ex_f = static_cast<Field2D *>( fields->filter_->Ex_[0] );
    Field2D *Ey_f = static_cast<Field2D *>( fields->filter_->Ey_[0] );
    Field2D *Ez_f = static_cast<Field2D *>( fields->filter_->Ez_[0] );
    Field2D *Ex_m1 = static_cast<Field2D *>( fields->filter_->Ex_[1] );
    Field2D *Ey_m1 = static_cast<Field2D *>( fields->filter_->Ey_[1] );
    Field2D *Ez_m1 = static_cast<Field2D *>( fields->filter_->Ez_[1] );
    Field2D *Ex_m2 = static_cast<Field2D *>( fields->filter_->Ex_[2] );
    Field2D *Ey_m2 = static_cast<Field2D *>( fields->filter_->Ey_[2] );
    Field2D *Ez_m2 = static_cast<Field2D *>( fields->filter_->Ez_[2] );
    
    double adv = 0.;
    
    // Electric field Ex^(d,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
        
            adv             = -dt*( *Jx2D )( i, j ) + dt_ov_dy * ( ( *Bz2D )( i, j+1 ) - ( *Bz2D )( i, j ) );
            // advance electric field
            ( *Ex2D )( i, j )   += adv;
            // compute the time-filtered field
            ( *Ex_f )( i, j )  = alpha*( *Ex2D )( i, j ) + beta*adv + delta*( ( *Ex_m1 )( i, j )+ftheta*( *Ex_m2 )( i, j ) );
            // update Ex_m2 and Ex_m1
            ( *Ex_m2 )( i, j )   = ( *Ex_m1 )( i, j ) - ftheta*( *Ex_m2 )( i, j );
            ( *Ex_m1 )( i, j )   = ( *Ex2D )( i, j )  - adv;
            
        }
    }
    
    
    // Electric field Ey^(p,d)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
        
            adv             = -dt*( *Jy2D )( i, j ) - dt_ov_dx * ( ( *Bz2D )( i+1, j ) - ( *Bz2D )( i, j ) );
            // advance electric field
            ( *Ey2D )( i, j )   += adv;
            // compute the time-filtered field
            ( *Ey_f )( i, j )  = alpha*( *Ey2D )( i, j ) + beta*adv + delta*( ( *Ey_m1 )( i, j )+ftheta*( *Ey_m2 )( i, j ) );
            // update Ex_m2 and Ex_m1
            ( *Ey_m2 )( i, j )   = ( *Ey_m1 )( i, j ) - ftheta*( *Ey_m2 )( i, j ) ;
            ( *Ey_m1 )( i, j )   = ( *Ey2D )( i, j )  - adv;
            
        }
    }
    
    
    // Electric field Ez^(p,p)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
        
            adv             = -dt*( *Jz2D )( i, j )
                              +                 dt_ov_dx * ( ( *By2D )( i+1, j ) - ( *By2D )( i, j ) )
                              -                 dt_ov_dy * ( ( *Bx2D )( i, j+1 ) - ( *Bx2D )( i, j ) );
            // advance electric field
            ( *Ez2D )( i, j )   += adv;
            // compute the time-filtered field
            ( *Ez_f )( i, j )  = alpha*( *Ez2D )( i, j ) + beta*adv + delta*( ( *Ez_m1 )( i, j )+ftheta*( *Ez_m2 )( i, j ) );
            // update Ex_m2 and Ex_m1
            ( *Ez_m2 )( i, j )   = ( *Ez_m1 )( i, j ) - ftheta * ( *Ez_m2 )( i, j ) ;
            ( *Ez_m1 )( i, j )   = ( *Ez2D )( i, j ) - adv ;
            
        }
    }
    
}

