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

ElectroMagnBC1D_SM::ElectroMagnBC1D_SM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC1D( params, patch, _min_max )
{
    // Parameters for the Silver-Mueller boundary conditions
    Alpha_SM = 2./( 1.+dt_ov_dx );
    Beta_SM  = ( dt_ov_dx-1. )/( 1.+dt_ov_dx );
    Gamma_SM = 4./( 1.+dt_ov_dx );
    
    By_val = 0.;
    Bz_val = 0.;
    
}

ElectroMagnBC1D_SM::~ElectroMagnBC1D_SM()
{
}

void ElectroMagnBC1D_SM::save_fields( Field *my_field, Patch *patch )
{
    Field1D *field1D=static_cast<Field1D *>( my_field );
    // Bx^(p) is not saved as it is defined on the primal grid and thus can be computed
    // we save only the field By and Bz that are computed on the dual grid
    
    if( min_max == 0 && patch->isXmin() ) {
        if( field1D->name=="By" ) {
            By_val = ( *my_field )( 0 );
        } else if( field1D->name=="Bz" ) {
            Bz_val = ( *my_field )( 0 );
        }
    } else if( min_max == 1 && patch->isXmax() ) {
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
    if( min_max == 0 && patch->isXmin() ) {
    
        //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D *Ey1D   = static_cast<Field1D *>( EMfields->Ey_ );
        Field1D *Ez1D   = static_cast<Field1D *>( EMfields->Ez_ );
        Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
        Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
        
        // Lasers
        double byL=0.0, bzL=0.0;
        vector<double> pos( 1 );
        pos[0] = 0.;
        for( unsigned int ilaser=0; ilaser<vecLaser.size(); ilaser++ ) {
            byL += vecLaser[ilaser]->getAmplitude0( pos, time_dual, 0, 0 );
            bzL += vecLaser[ilaser]->getAmplitude1( pos, time_dual, 0, 0 );
        }
        
        // Apply Silver-Mueller EM boundary condition at x=xmin
        ( *By1D )( 0 ) =  Alpha_SM*( *Ez1D )( 0 ) + Beta_SM*( ( *By1D )( 1 )-By_val ) + Gamma_SM*byL+By_val;
        ( *Bz1D )( 0 ) = -Alpha_SM*( *Ey1D )( 0 ) + Beta_SM*( ( *Bz1D )( 1 )-Bz_val ) + Gamma_SM*bzL+Bz_val;
        
    } else if( min_max == 1 && patch->isXmax() ) {
        //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D *Ey1D   = static_cast<Field1D *>( EMfields->Ey_ );
        Field1D *Ez1D   = static_cast<Field1D *>( EMfields->Ez_ );
        Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
        Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
        
        // Lasers
        double byR=0.0, bzR=0.0;
        vector<double> pos( 1 );
        pos[0] = 0.;
        for( unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++ ) {
            byR += vecLaser[ilaser]->getAmplitude0( pos, time_dual, 0, 0 );
            bzR += vecLaser[ilaser]->getAmplitude1( pos, time_dual, 0, 0 );
        }
        
        // Silver-Mueller boundary conditions (right)
        ( *By1D )( nx_d-1 ) = -Alpha_SM*( *Ez1D )( nx_p-1 )+ Beta_SM*( ( *By1D )( nx_d-2 )-By_val ) + Gamma_SM*byR+By_val;
        ( *Bz1D )( nx_d-1 ) =  Alpha_SM*( *Ey1D )( nx_p-1 )+ Beta_SM*( ( *Bz1D )( nx_d-2 )-Bz_val ) + Gamma_SM*bzR+Bz_val;
    }
    
}
