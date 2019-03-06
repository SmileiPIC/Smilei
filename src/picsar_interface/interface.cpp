#include "interface.h"


void copy_field_3d( Field3D *out, Field3D *in )
{
    unsigned int n1, n2, n3;
    unsigned int i, j, k;
    n1 = ( ( in->dims_[0] ) < ( out->dims_[0] ) ? ( in->dims_[0] ) : ( out->dims_[0] ) );
    n2 = ( ( in->dims_[1] ) < ( out->dims_[1] ) ? ( in->dims_[1] ) : ( out->dims_[1] ) );
    n3 = ( ( in->dims_[2] ) < ( out->dims_[2] ) ? ( in->dims_[2] ) : ( out->dims_[2] ) );
    #pragma omp parallel for collapse(3) private(i ,j ,k) schedule(runtime)
    for( i=0; i<n1; i++ ) {
        for( j=0; j<n2; j++ ) {
            for( k=0; k<n3; k++ ) {
                ( *out )( i, j, k ) = ( *in )( i, j, k );
            }
        }
    }
}

void copy_field_2d( Field2D *out, Field2D *in )
{
    unsigned int n1, n2;
    unsigned int i, j;
    n1 = ( ( in->dims_[0] ) < ( out->dims_[0] ) ? ( in->dims_[0] ) : ( out->dims_[0] ) );
    n2 = ( ( in->dims_[1] ) < ( out->dims_[1] ) ? ( in->dims_[1] ) : ( out->dims_[1] ) );
    #pragma omp parallel for collapse(2) private(i , j )  schedule(runtime)
    for( i=0; i<n1; i++ ) {
        for( j=0; j<n2; j++ ) {
            ( *out )( i, j ) = ( *in )( i, j );
        }
    }
}


void duplicate_field_into_pxr( ElectroMagn *fields )
{
    int nDim_field = fields->Ex_->dims_.size();
    if( nDim_field == 3 ) {
        Field3D *Ex3D_pxr = static_cast<Field3D *>( fields->Ex_pxr );
        Field3D *Ey3D_pxr = static_cast<Field3D *>( fields->Ey_pxr );
        Field3D *Ez3D_pxr = static_cast<Field3D *>( fields->Ez_pxr );
        Field3D *Bx3D_pxr = static_cast<Field3D *>( fields->Bx_pxr );
        Field3D *By3D_pxr = static_cast<Field3D *>( fields->By_pxr );
        Field3D *Bz3D_pxr = static_cast<Field3D *>( fields->Bz_pxr );
        Field3D *Jx3D_pxr = static_cast<Field3D *>( fields->Jx_pxr );
        Field3D *Jy3D_pxr = static_cast<Field3D *>( fields->Jy_pxr );
        Field3D *Jz3D_pxr = static_cast<Field3D *>( fields->Jz_pxr );
        Field3D *rho3D_pxr = static_cast<Field3D *>( fields->rho_pxr );
        Field3D *rhoold3D_pxr = static_cast<Field3D *>( fields->rhoold_pxr );
        
        Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
        Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
        Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
        Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
        Field3D *By3D = static_cast<Field3D *>( fields->By_ );
        Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
        Field3D *Jx3D = static_cast<Field3D *>( fields->Jx_ );
        Field3D *Jy3D = static_cast<Field3D *>( fields->Jy_ );
        Field3D *Jz3D = static_cast<Field3D *>( fields->Jz_ );
        Field3D *rho3D = static_cast<Field3D *>( fields->rho_ );
        Field3D *rhoold3D = static_cast<Field3D *>( fields->rhoold_ );
        copy_field_3d( Ex3D_pxr, Ex3D );
        copy_field_3d( Ey3D_pxr, Ey3D );
        copy_field_3d( Ez3D_pxr, Ez3D );
        copy_field_3d( Bx3D_pxr, Bx3D );
        copy_field_3d( By3D_pxr, By3D );
        copy_field_3d( Bz3D_pxr, Bz3D );
        copy_field_3d( Jx3D_pxr, Jx3D );
        copy_field_3d( Jy3D_pxr, Jy3D );
        copy_field_3d( Jz3D_pxr, Jz3D );
        copy_field_3d( rho3D_pxr, rho3D );
        copy_field_3d( rhoold3D_pxr, rhoold3D );
    }
    if( nDim_field ==2 ) {
        Field2D *Ex2D_pxr = static_cast<Field2D *>( fields->Ex_pxr );
        Field2D *Ey2D_pxr = static_cast<Field2D *>( fields->Ey_pxr );
        Field2D *Ez2D_pxr = static_cast<Field2D *>( fields->Ez_pxr );
        Field2D *Bx2D_pxr = static_cast<Field2D *>( fields->Bx_pxr );
        Field2D *By2D_pxr = static_cast<Field2D *>( fields->By_pxr );
        Field2D *Bz2D_pxr = static_cast<Field2D *>( fields->Bz_pxr );
        Field2D *Jx2D_pxr = static_cast<Field2D *>( fields->Jx_pxr );
        Field2D *Jy2D_pxr = static_cast<Field2D *>( fields->Jy_pxr );
        Field2D *Jz2D_pxr = static_cast<Field2D *>( fields->Jz_pxr );
        Field2D *rho2D_pxr = static_cast<Field2D *>( fields->rho_pxr );
        Field2D *rhoold2D_pxr = static_cast<Field2D *>( fields->rhoold_pxr );
        
        Field2D *Ex2D = static_cast<Field2D *>( fields->Ex_ );
        Field2D *Ey2D = static_cast<Field2D *>( fields->Ey_ );
        Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( fields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
        Field2D *Jx2D = static_cast<Field2D *>( fields->Jx_ );
        Field2D *Jy2D = static_cast<Field2D *>( fields->Jy_ );
        Field2D *Jz2D = static_cast<Field2D *>( fields->Jz_ );
        Field2D *rho2D = static_cast<Field2D *>( fields->rho_ );
        Field2D *rhoold2D = static_cast<Field2D *>( fields->rhoold_ );
        copy_field_2d( Ex2D_pxr, Ex2D );
        copy_field_2d( Ey2D_pxr, Ey2D );
        copy_field_2d( Ez2D_pxr, Ez2D );
        copy_field_2d( Bx2D_pxr, Bx2D );
        copy_field_2d( By2D_pxr, By2D );
        copy_field_2d( Bz2D_pxr, Bz2D );
        copy_field_2d( Jx2D_pxr, Jx2D );
        copy_field_2d( Jy2D_pxr, Jy2D );
        copy_field_2d( Jz2D_pxr, Jz2D );
        copy_field_2d( rho2D_pxr, rho2D );
        copy_field_2d( rhoold2D_pxr, rhoold2D );
    }
}

void duplicate_field_into_smilei( ElectroMagn *fields )
{
    int nDim_field = fields->Ex_->dims_.size();
    if( nDim_field ==3 ) {
        Field3D *Ex3D_pxr = static_cast<Field3D *>( fields->Ex_pxr );
        Field3D *Ey3D_pxr = static_cast<Field3D *>( fields->Ey_pxr );
        Field3D *Ez3D_pxr = static_cast<Field3D *>( fields->Ez_pxr );
        Field3D *Bx3D_pxr = static_cast<Field3D *>( fields->Bx_pxr );
        Field3D *By3D_pxr = static_cast<Field3D *>( fields->By_pxr );
        Field3D *Bz3D_pxr = static_cast<Field3D *>( fields->Bz_pxr );
        
        
        Field3D *Ex3D = static_cast<Field3D *>( fields->Ex_ );
        Field3D *Ey3D = static_cast<Field3D *>( fields->Ey_ );
        Field3D *Ez3D = static_cast<Field3D *>( fields->Ez_ );
        Field3D *Bx3D = static_cast<Field3D *>( fields->Bx_ );
        Field3D *By3D = static_cast<Field3D *>( fields->By_ );
        Field3D *Bz3D = static_cast<Field3D *>( fields->Bz_ );
        
        
        copy_field_3d( Ex3D, Ex3D_pxr );
        copy_field_3d( Ey3D, Ey3D_pxr );
        copy_field_3d( Ez3D, Ez3D_pxr );
        copy_field_3d( Bx3D, Bx3D_pxr );
        copy_field_3d( By3D, By3D_pxr );
        copy_field_3d( Bz3D, Bz3D_pxr );
        
    }
    if( nDim_field ==2 ) {
        Field2D *Ex2D_pxr = static_cast<Field2D *>( fields->Ex_pxr );
        Field2D *Ey2D_pxr = static_cast<Field2D *>( fields->Ey_pxr );
        Field2D *Ez2D_pxr = static_cast<Field2D *>( fields->Ez_pxr );
        Field2D *Bx2D_pxr = static_cast<Field2D *>( fields->Bx_pxr );
        Field2D *By2D_pxr = static_cast<Field2D *>( fields->By_pxr );
        Field2D *Bz2D_pxr = static_cast<Field2D *>( fields->Bz_pxr );
        
        
        Field2D *Ex2D = static_cast<Field2D *>( fields->Ex_ );
        Field2D *Ey2D = static_cast<Field2D *>( fields->Ey_ );
        Field2D *Ez2D = static_cast<Field2D *>( fields->Ez_ );
        Field2D *Bx2D = static_cast<Field2D *>( fields->Bx_ );
        Field2D *By2D = static_cast<Field2D *>( fields->By_ );
        Field2D *Bz2D = static_cast<Field2D *>( fields->Bz_ );
        
        
        
        copy_field_2d( Ex2D, Ex2D_pxr );
        copy_field_2d( Ey2D, Ey2D_pxr );
        copy_field_2d( Ez2D, Ez2D_pxr );
        copy_field_2d( Bx2D, Bx2D_pxr );
        copy_field_2d( By2D, By2D_pxr );
        copy_field_2d( Bz2D, Bz2D_pxr );
        
    }
}

