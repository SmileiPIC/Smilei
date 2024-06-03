
#include "MF_Solver3D_Yee.h"

#include "ElectroMagn.h"
#include "Field3D.h"

MF_Solver3D_Yee::MF_Solver3D_Yee( Params &params )
    : Solver3D( params )
{
    isEFilterApplied = params.Friedman_filter;
}

MF_Solver3D_Yee::~MF_Solver3D_Yee()
{
    // EMPTY
}

void MF_Solver3D_Yee::operator()( ElectroMagn *fields )
{
    // Static-cast of the fields
    double *const __restrict__ Bx3D       = fields->Bx_->data();
    double *const __restrict__ By3D       = fields->By_->data();
    double *const __restrict__ Bz3D       = fields->Bz_->data();

    const unsigned int nx_p = fields->dimPrim[0];
    const unsigned int nx_d = fields->dimDual[0];
    const unsigned int ny_p = fields->dimPrim[1];
    const unsigned int ny_d = fields->dimDual[1];
    const unsigned int nz_p = fields->dimPrim[2];
    const unsigned int nz_d = fields->dimDual[2];
    //double *__restrict__ Ex3D ;
    const double * __restrict__ Ex3D = isEFilterApplied ? fields->filter_->Ex_[0]->data() : fields->Ex_->data();
    const double * __restrict__ Ey3D = isEFilterApplied ? fields->filter_->Ey_[0]->data() : fields->Ey_->data();
    const double * __restrict__ Ez3D = isEFilterApplied ? fields->filter_->Ez_[0]->data() : fields->Ez_->data();

    // Magnetic field Bx^(p,d,d)
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    const int sizeofEx = fields->Ex_->number_of_points_;
    const int sizeofEy = fields->Ey_->number_of_points_;
    const int sizeofEz = fields->Ez_->number_of_points_;
    const int sizeofBx = fields->Bx_->number_of_points_;
    const int sizeofBy = fields->By_->number_of_points_;
    const int sizeofBz = fields->Bz_->number_of_points_;

    #pragma acc parallel present( Bx3D[0:sizeofBx], Ey3D[0:sizeofEy], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=0 ; i<nx_p;  i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                Bx3D[ i*(ny_d*nz_d) + j*(nz_d) + k ] += -dt_ov_dy * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ i*(ny_p*nz_d) + (j-1)*(nz_d) + k   ] )
                                                     +   dt_ov_dz * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ i*(ny_d*nz_p) +  j   *(nz_p) + k-1 ] );
            }
        }
    }

    // Magnetic field By^(d,p,d)
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( By3D[0:sizeofBy], Ex3D[0:sizeofEx], Ez3D[0:sizeofEz] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=0 ; j<ny_p ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=1 ; k<nz_d-1 ; k++ ) {
                By3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] += -dt_ov_dz * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + j*(nz_p) + k-1 ] )
                                                     +   dt_ov_dx * ( Ez3D[ i*(ny_p*nz_d) + j*(nz_d) + k ] - Ez3D[ (i-1)*(ny_p*nz_d) + j*(nz_d) + k   ] );
            }
        }
    }

    // Magnetic field Bz^(d,d,p)
#if defined( SMILEI_ACCELERATOR_GPU_OACC )
    #pragma acc parallel present( Bz3D[0:sizeofBz], Ex3D[0:sizeofEx], Ey3D[0:sizeofEy] )
    #pragma acc loop gang
#elif defined( SMILEI_ACCELERATOR_GPU_OMP )
    #pragma omp target
    #pragma omp teams distribute parallel for collapse( 3 )
#endif
    for( unsigned int i=1 ; i<nx_d-1 ; i++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
        #pragma acc loop worker
#endif
        for( unsigned int j=1 ; j<ny_d-1 ; j++ ) {
#ifdef SMILEI_ACCELERATOR_GPU_OACC
            #pragma acc loop vector
#endif
            for( unsigned int k=0 ; k<nz_p ; k++ ) {
                Bz3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] += -dt_ov_dx * ( Ey3D[ i*(ny_d*nz_p) + j*(nz_p) + k ] - Ey3D[ (i-1)*(ny_d*nz_p) +  j   *(nz_p) + k ] )
                                                     +   dt_ov_dy * ( Ex3D[ i*(ny_p*nz_p) + j*(nz_p) + k ] - Ex3D[  i   *(ny_p*nz_p) + (j-1)*(nz_p) + k ] );
            }
        }
    }
}
