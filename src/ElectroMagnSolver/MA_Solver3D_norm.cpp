
#include "MA_Solver3D_norm.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MA_Solver3D_norm::MA_Solver3D_norm( Params &params )
    : Solver3D( params )
{
}

MA_Solver3D_norm::~MA_Solver3D_norm()
{
}

void MA_Solver3D_norm::operator()( ElectroMagn *fields )
{

    // Static-cast of the fields
    double *Ex3D = &(fields->Ex_->data_[0]);
    double *Ey3D = &(fields->Ey_->data_[0]);
    double *Ez3D = &(fields->Ez_->data_[0]);
    double *Bx3D = &(fields->Bx_->data_[0]);
    double *By3D = &(fields->By_->data_[0]);
    double *Bz3D = &(fields->Bz_->data_[0]);
    double *Jx3D = &(fields->Jx_->data_[0]);
    double *Jy3D = &(fields->Jy_->data_[0]);
    double *Jz3D = &(fields->Jz_->data_[0]);
    
    // Electric field Ex^(d,p,p)
    for( unsigned int i=0 ; i<nx_d ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] += -dt*Jx3D[ i*(ny_p*nz_p) + j*(nz_p) + k ]
                    +                 dt_ov_dy * ( Bz3D[ i*(ny_d*nz_p) + (j+1)*(nz_p) + k   ] - Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] )
                    -                 dt_ov_dz * ( By3D[ i*(ny_p*nz_d) +  j   *(nz_d) + k+1 ] - By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] );
            }
        }
    }
    
    // Electric field Ey^(p,d,p)
    for( unsigned int i=0 ; i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_d ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt*Jy3D[ i*(ny_d*nz_p) + j*(nz_p) + k ]
                    -                  dt_ov_dx * ( Bz3D[ (i+1)*(ny_d*nz_p) + j*(nz_p) + k   ] - Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] )
                    +                  dt_ov_dz * ( Bx3D[  i   *(ny_d*nz_d) + j*(nz_d) + k+1 ] - Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] );
            }
        }
    }
    
    // Electric field Ez^(p,p,d)
    for( unsigned int i=0 ;  i<nx_p ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=0 ; k<nz_d ; k++ ) {
                Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt*Jz3D[ i*(ny_p*nz_d) + j*(nz_d) + k ]
                    +                  dt_ov_dx * ( By3D[ (i+1)*(ny_p*nz_d) +  j   *(nz_d) + k ] - By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] )
                    -                  dt_ov_dy * ( Bx3D[  i   *(ny_d*nz_d) + (j+1)*(nz_d) + k ] - Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] );
            }
        }
    }
    
}

