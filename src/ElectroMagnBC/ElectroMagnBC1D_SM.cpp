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
    Alpha_ = 2. / ( 1. + dt_ov_d[0] );
    Beta_  = ( dt_ov_d[0] - 1. ) / ( 1. + dt_ov_d[0] );
    Gamma_ = 4. / ( 1. + dt_ov_d[0] );
    
    By_val_ = 0.;
    Bz_val_ = 0.;
    
    sign_ = (double) (i_boundary_ % 2) *2 - 1.; // -1 or 1 for min or max
    
    if( i_boundary == 0 ) {
        iE_ = 0;
        iB_ = 0;
        iB_old_ = 1;
    } else {
        iE_ = n_p[0] - 1;
        iB_ = n_d[0] - 1;
        iB_old_ = iB_ - 1;
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
            By_val_ = ( *my_field )( 0 );
        } else if( field1D->name=="Bz" ) {
            Bz_val_ = ( *my_field )( 0 );
        }
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
        if( field1D->name=="By" ) {
            By_val_ = ( *my_field )( field1D->dims()[0]-1 );
        } else if( field1D->name=="Bz" ) {
            Bz_val_ = ( *my_field )( field1D->dims()[0]-1 );
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
        /*Field1D *Ey1D   = static_cast<Field1D *>( EMfields->Ey_ );
        Field1D *Ez1D   = static_cast<Field1D *>( EMfields->Ez_ );
        Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
        Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );*/
        
        const Field  *E[3]{ EMfields->Ex_, EMfields->Ey_, EMfields->Ez_ };
        const Field  *B[3]{ EMfields->Bx_, EMfields->By_, EMfields->Bz_ };
        const double *const __restrict__ E1 = E[1]->data_;
        const double *const __restrict__ E2 = E[2]->data_;
        double *const __restrict__ B1       = B[1]->data_;
        double *const __restrict__ B2       = B[2]->data_;
        // Lasers
        double by = 0., bz = 0.;
        vector<double> pos( 1 );
        pos[0] = 0.;
        for( unsigned int ilaser=0; ilaser<vecLaser.size(); ilaser++ ) {
            by += vecLaser[ilaser]->getAmplitude0( pos, time_dual, 0, 0 );
            bz += vecLaser[ilaser]->getAmplitude1( pos, time_dual, 0, 0 );
        }
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        const int sizeofE1 = E[1]->number_of_points_;
        const int sizeofE2 = E[2]->number_of_points_;
        const int sizeofB1 = B[1]->number_of_points_;
        const int sizeofB2 = B[2]->number_of_points_;
#endif
        // Apply Silver-Mueller EM boundary condition at x=xmin or xmax
        
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc parallel present(E1[0:sizeofE1],E2[0:sizeofE2],B1[0:sizeofB1],B2[0:sizeofB2])
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
        #pragma omp target
	#pragma omp teams distribute parallel for
	for (int i=0; i<1;++i)
#endif
        {
            //( *By1D )( iB_ ) = -sign_*Alpha_*( *Ez1D )( iE_ ) + Beta_*( ( *By1D )( iB_old_ )-By_val_ ) + Gamma_*by + By_val_;
            //( *Bz1D )( iB_ ) =  sign_*Alpha_*( *Ey1D )( iE_ ) + Beta_*( ( *Bz1D )( iB_old_ )-Bz_val_ ) + Gamma_*bz + Bz_val_;
            B1[ iB_ ] = -sign_ * Alpha_ * E2[iE_] + Beta_ * ( B1[iB_old_] - By_val_) + Gamma_ * by + By_val_;
            B2[ iB_ ] =  sign_ * Alpha_ * E1[iE_] + Beta_ * ( B2[iB_old_] - Bz_val_) + Gamma_ * bz + Bz_val_;
        }
    }
    
}
