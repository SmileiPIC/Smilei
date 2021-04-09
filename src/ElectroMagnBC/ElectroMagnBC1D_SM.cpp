#include "ElectroMagnBC1D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC1D_SM::ElectroMagnBC1D_SM( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC1D( params, patch, i_boundary )
{
    // Parameters for the Silver-Mueller boundary conditions
    Alpha = 2./( 1.+dt_ov_d[0] );
    Beta  = ( dt_ov_d[0]-1. )/( 1.+dt_ov_d[0] );
    Gamma = 4./( 1.+dt_ov_d[0] );
    
    By_val = 0.;
    Bz_val = 0.;
    
    sign_ = (double) (i_boundary_ % 2) *2 - 1.; // -1 or 1 for min or max
    
    if( i_boundary == 0 ) {
        iE = 0;
        iB = 0;
        iB_old = 1;
    } else {
        iE = n_p[0] - 1;
        iB = n_d[0] - 1;
        iB_old = iB - 1;
    }
    
}

ElectroMagnBC1D_SM::~ElectroMagnBC1D_SM()
{
}

void ElectroMagnBC1D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field1D *field1D=static_cast<Field1D *>( my_field );
    // Bx^(p) is not saved as it is defined on the primal grid and thus can be computed
    // we save only the field By and Bz that are computed on the dual grid
    
    if( i_boundary_ == 0 && patch->isXmin() ) {
        if( field1D->name=="By" ) {
            By_val = ( *my_field )( 0 );
        } else if( field1D->name=="Bz" ) {
            Bz_val = ( *my_field )( 0 );
        }
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
        if( field1D->name=="By" ) {
            By_val = ( *my_field )( field1D->dims()[0]-1 );
        } else if( field1D->name=="Bz" ) {
            Bz_val = ( *my_field )( field1D->dims()[0]-1 );
        }
        
    }
}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( patch->isBoundary( i_boundary_ ) ) {
    
        //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D *Ey1D   = static_cast<Field1D *>( EMfields->Ey_ );
        Field1D *Ez1D   = static_cast<Field1D *>( EMfields->Ez_ );
        Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
        Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
        
        // Lasers
        double by = 0., bz = 0.;
        vector<double> pos( 1 );
        pos[0] = 0.;
        for( unsigned int ilaser=0; ilaser<vecLaser.size(); ilaser++ ) {
            by += vecLaser[ilaser]->getAmplitude0( pos, time_dual, 0, 0 );
            bz += vecLaser[ilaser]->getAmplitude1( pos, time_dual, 0, 0 );
        }
        
        // Apply Silver-Mueller EM boundary condition at x=xmin or xmax
        
        ( *By1D )( iB ) = -sign_*Alpha*( *Ez1D )( iE ) + Beta*( ( *By1D )( iB_old )-By_val ) + Gamma*by + By_val;
        ( *Bz1D )( iB ) =  sign_*Alpha*( *Ey1D )( iE ) + Beta*( ( *Bz1D )( iB_old )-Bz_val ) + Gamma*bz + Bz_val;
        
    }
    
}
