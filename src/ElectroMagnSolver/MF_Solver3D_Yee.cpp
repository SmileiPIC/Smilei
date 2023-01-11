
#include "MF_Solver3D_Yee.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MF_Solver3D_Yee::MF_Solver3D_Yee( Params &params )
    : Solver3D( params )
{
}

MF_Solver3D_Yee::~MF_Solver3D_Yee()
{
}

void MF_Solver3D_Yee::operator()( ElectroMagn *fields )
{
    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    const unsigned int nz_p = fields->dimPrim[2];
    const unsigned int nz_d = fields->dimDual[2];
    // Static-cast of the fields
    double *Ex3D = &(fields->Ex_->data_[0]);
    double *Ey3D = &(fields->Ey_->data_[0]);
    double *Ez3D = &(fields->Ez_->data_[0]);
    double *Bx3D = &(fields->Bx_->data_[0]);
    double *By3D = &(fields->By_->data_[0]);
    double *Bz3D = &(fields->Bz_->data_[0]);
    
    // Magnetic field Bx^(p,d,d)
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
            }
        }
    }
    
}

