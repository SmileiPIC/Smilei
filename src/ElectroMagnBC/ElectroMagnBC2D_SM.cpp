#include "ElectroMagnBC2D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC2D_SM::ElectroMagnBC2D_SM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC2D( params, patch, i_boundary )
{
    
    // Calculate axes indices
    axis0_ = i_boundary_ / 2; // axis normal to border
    axis1_ = axis0_ == 0 ? 1 : 0; // other axis (tangent to border)
    sign_ = (double) (i_boundary_ % 2) *2 - 1.; // -1 or 1 for min or max
    
    // Index where to set the field along axis 0
    iB_.resize( 3 );
    if( sign_ < 0 ) {
        iB_[0] = 0;
        iB_[1] = 0;
        iB_[2] = 0;
    } else {
        iB_[axis0_] = n_p[axis0_] - 1;
        iB_[axis1_] = n_d[axis0_] - 1;
        iB_[2     ] = n_d[axis0_] - 1;
    }
    
    // Buffers to save B field
    B_val.resize( 3 );
    if( patch->isBoundary( i_boundary_ ) ) {
        B_val[axis0_].resize( n_d[axis1_], 0. ); // dual in the first direction
        B_val[axis1_].resize( n_p[axis1_], 0. ); // primal in the other direction
        B_val[2     ].resize( n_d[axis1_], 0. ); // dual in the other direction
    }
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------
    
    vector<double> K = params.EM_BCs_k[i_boundary_];
    double Knorm = sqrt( K[0]*K[0] + K[1]*K[1] ) ;
    double omega = 1.;
    double k0 = omega*K[axis0_] / Knorm;
    double k1 = omega*K[axis1_] / Knorm;
    
    double factor = 1.0 / ( k0 - sign_ * dt_ov_d[axis0_] );
    Alpha_   = 2.0 * factor;
    Beta_    = - ( k0 + sign_ * dt_ov_d[axis0_] ) * factor;
    Gamma_   = 4.0 * k0 * factor;
    Delta_   = - ( k1 + dt_ov_d[axis1_] ) * factor;
    Epsilon_ = - ( k1 - dt_ov_d[axis1_] ) * factor;
}


void ElectroMagnBC2D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field2D *field2D=static_cast<Field2D *>( my_field );
    
    if( patch->isBoundary( i_boundary_ ) ) {
        
        unsigned int xyz = 0;
        if( field2D->name=="Bx" ) {
            xyz = 0;
        } else if( field2D->name=="By" ) {
            xyz = 1;
        } else if( field2D->name=="Bz" ) {
            xyz = 2;
        }
        
        if( axis0_ == 0 ) {
            for( unsigned int j=0; j<B_val[xyz].size(); j++ ) {
                B_val[xyz][j]=( *field2D )( iB_[xyz], j );
            }
        } else {
            for( unsigned int i=0; i<B_val[xyz].size(); i++ ) {
                B_val[xyz][i]=( *field2D )( i, iB_[xyz] );
            }
        }
        
    }
}

void ElectroMagnBC2D_SM::disableExternalFields()
{
    B_val[0].resize( 0 );
    B_val[1].resize( 0 );
    B_val[2].resize( 0 );
}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC2D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    
    if( patch->isBoundary( i_boundary_ ) ) {
        
        // Static cast of the fields
        vector<Field2D*> E( 3 );
        E[0] = static_cast<Field2D *>( EMfields->Ex_ );
        E[1] = static_cast<Field2D *>( EMfields->Ey_ );
        E[2] = static_cast<Field2D *>( EMfields->Ez_ );
        vector<Field2D*> B( 3 );
        B[0] = static_cast<Field2D *>( EMfields->Bx_ );
        B[1] = static_cast<Field2D *>( EMfields->By_ );
        B[2] = static_cast<Field2D *>( EMfields->Bz_ );
        
        // Lasers polarized along axis 1
        vector<double> b1( n_p[axis1_], 0. );
        vector<double> pos( 1 );
        if( ! vecLaser.empty() ) {
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_p[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    b1[j] += vecLaser[ilaser]->getAmplitude0( pos, time_dual, j, 0 );
                }
            }
        }
        if( axis0_ == 0 ) { // for By^(d,p)
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_p[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                ( *B[1] )( iB_[1], j )
                    = Alpha_  *   ( *E[2] )( iB_[0]      , j )
                    + Beta_   * ( ( *B[1] )( iB_[1]-sign_, j )-B_val[1][j] )
                    + Gamma_  * b1[j]
                    + Delta_  * ( ( *B[0] )( iB_[0], j+1 )-B_val[0][j+1] )
                    + Epsilon_* ( ( *B[0] )( iB_[0], j   )-B_val[0][j  ] )
                    + B_val[1][j];
            }
        } else { // for Bx^(p,d)
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_p[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                ( *B[0] )( j, iB_[0] )
                    = -Alpha_ *   ( *E[2] )( j, iB_[1]        )
                    + Beta_   * ( ( *B[0] )( j, iB_[0]-sign_ )-B_val[0][j] )
                    + Gamma_  * b1[j]
                    + Delta_  * ( ( *B[1] )( j+1, iB_[1] )-B_val[1][j+1] )
                    + Epsilon_* ( ( *B[1] )( j  , iB_[1] )-B_val[1][j  ] )
                    + B_val[0][j];
            }
        }
        
        // Lasers polarized along axis 2
        vector<double> b2( n_d[axis1_], 0. );
        if( ! vecLaser.empty() ) {
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_d[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                pos[0] = patch->getDomainLocalMin( axis1_ ) + ( ( int )j - 0.5 - ( int )EMfields->oversize[axis1_] )*d[axis1_];
                for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
                    b2[j] += vecLaser[ilaser]->getAmplitude1( pos, time_dual, j, 0 );
                }
            }
        }
        // for Bz^(d,d)
        if( axis0_ == 0 ) {
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_d[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                ( *B[2] )( iB_[2], j )
                    = -Alpha_ *   ( *E[1] )( iB_[0]      , j )
                    + Beta_   * ( ( *B[2] )( iB_[2]-sign_, j )- B_val[2][j] )
                    + Gamma_  * b2[j]
                    + B_val[2][j];
            }
        } else {
            for( unsigned int j=patch->isBoundary(axis1_,0) ; j<n_d[axis1_]-patch->isBoundary(axis1_,1) ; j++ ) {
                ( *B[2] )( j, iB_[2] )
                    = Alpha_ *   ( *E[0] )( j, iB_[1]       )
                    + Beta_  * ( ( *B[2] )( j, iB_[2]-sign_ )- B_val[2][j] )
                    + Gamma_ * b2[j]
                    + B_val[2][j];
            }
        }
        
    }
}
