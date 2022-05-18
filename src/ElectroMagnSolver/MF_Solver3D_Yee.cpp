
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
    // Static-cast of the fields
    const double *const __restrict__ Ex3D = &( fields->Ex_->data_[0] );
    const double *const __restrict__ Ey3D = &( fields->Ey_->data_[0] );
    const double *const __restrict__ Ez3D = &( fields->Ez_->data_[0] );
    double *const __restrict__ Bx3D       = &( fields->Bx_->data_[0] );
    double *const __restrict__ By3D       = &( fields->By_->data_[0] );
    double *const __restrict__ Bz3D       = &( fields->Bz_->data_[0] );

    const int sizeofEx = fields->Ex_->globalDims_;
    const int sizeofEy = fields->Ey_->globalDims_;
    const int sizeofEz = fields->Ez_->globalDims_;
    const int sizeofBx = fields->Bx_->globalDims_;
    const int sizeofBy = fields->By_->globalDims_;
    const int sizeofBz = fields->Bz_->globalDims_;

    // Magnetic field Bx^(p,d,d)
#if defined( _GPU )
    #pragma acc parallel present( Bx3D[0:sizeofBx], Ey3D[0:sizeofEy], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
#ifdef _GPU
        #pragma acc loop worker
#endif
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
            }
        }
    }
    
    // Magnetic field By^(d,p,d)
#if defined( _GPU )
    #pragma acc parallel present( By3D[0:sizeofBy], Ex3D[0:sizeofEx], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
#ifdef _GPU
        #pragma acc loop worker
#endif
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
            }
        }
    }
    
    // Magnetic field Bz^(d,d,p)
#if defined( _GPU )
    #pragma acc parallel present( Bz3D[0:sizeofBz], Ex3D[0:sizeofEx], Ey3D[0:sizeofEy] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams /* num_teams(xxx) thread_limit(xxx) */ // TODO(Etienne M): WG/WF tuning
    #pragma omp distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
#ifdef _GPU
        #pragma acc loop worker
#endif
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
#ifdef _GPU
            #pragma acc loop vector
#endif
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
            }
        }
    }
    
}

