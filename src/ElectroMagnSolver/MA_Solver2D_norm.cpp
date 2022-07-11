
#include "MA_Solver2D_norm.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MA_Solver2D_norm::MA_Solver2D_norm( Params &params )
    : Solver2D( params )
{
    // EMPTY
}

MA_Solver2D_norm::~MA_Solver2D_norm()
{
    // EMPTY
}

void MA_Solver2D_norm::operator()( ElectroMagn *fields )
{
    double *const __restrict__ Ex2D       = fields->Ex_->data(); // [x * ny_p + y] : dual in x   primal in y,z
    double *const __restrict__ Ey2D       = fields->Ey_->data(); // [x * ny_d + y] : dual in y   primal in x,z
    double *const __restrict__ Ez2D       = fields->Ez_->data(); // [x * ny_p + y] : dual in z   primal in x,y
    const double *const __restrict__ Bx2D = fields->Bx_->data(); // [x * ny_d + y] : dual in y,z primal in x
    const double *const __restrict__ By2D = fields->By_->data(); // [x * ny_p + y] : dual in x,z primal in y
    const double *const __restrict__ Bz2D = fields->Bz_->data(); // [x * ny_d + y] : dual in x,y primal in z
    const double *const __restrict__ Jx2D = fields->Jx_->data(); // [x * ny_p + y] : dual in x   primal in y,z
    const double *const __restrict__ Jy2D = fields->Jy_->data(); // [x * ny_d + y] : dual in y   primal in x,z
    const double *const __restrict__ Jz2D = fields->Jz_->data(); // [x * ny_p + y] : dual in z   primal in x,y

    // Electric field Ex^(d,p)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    #pragma omp target map( tofrom                                  \
                            : Ex2D [0:( nx_d - 1 ) * ny_p + ny_p] ) \
        map( to                                                     \
             : Jx2D [0:( nx_d - 1 ) * ny_p + ny_p],                 \
               Bz2D [0:( nx_d - 1 ) * ny_d + ny_p + 1] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 0; x < nx_d; ++x ) {
        for( unsigned int y = 0; y < ny_p; ++y ) {
            Ex2D[x * ny_p + y] += -dt * Jx2D[x * ny_p + y] + dt_ov_dy * ( Bz2D[x * ny_d + y + 1] - Bz2D[x * ny_d + y] );
        }
    }

    // Electric field Ey^(p,d)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    #pragma omp target map( tofrom                                  \
                            : Ey2D [0:( nx_p - 1 ) * ny_d + ny_d] ) \
        map( to                                                     \
             : Jy2D [0:( nx_p - 1 ) * ny_d + ny_d],                 \
               Bz2D [0:( nx_p - 1 + 1 ) * ny_d + ny_d] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 0; x < nx_p; ++x ) {
        for( unsigned int y = 0; y < ny_d; ++y ) {
            Ey2D[x * ny_d + y] += -dt * Jy2D[x * ny_d + y] - dt_ov_dx * ( Bz2D[( x + 1 ) * ny_d + y] - Bz2D[x * ny_d + y] );
        }
    }

    // Electric field Ez^(p,p)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    #pragma omp target map( tofrom                                  \
                            : Ez2D [0:( nx_p - 1 ) * ny_p + ny_p] ) \
        map( to                                                     \
             : Jz2D [0:( nx_p - 1 ) * ny_p + ny_p],                 \
               Bx2D [0:( nx_p - 1 ) * ny_d + ny_p + 1],             \
               By2D [0:( nx_p - 1 + 1 ) * ny_p + ny_p] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 0; x < nx_p; ++x ) {
        for( unsigned int y = 0; y < ny_p; ++y ) {
            Ez2D[x * ny_p + y] += -dt * Jz2D[x * ny_p + y] +
                                  dt_ov_dx * ( By2D[( x + 1 ) * ny_p + y] - By2D[x * ny_p + y] ) -
                                  dt_ov_dy * ( Bx2D[x * ny_d + y + 1] - Bx2D[x * ny_d + y] );
        }
    }
}
