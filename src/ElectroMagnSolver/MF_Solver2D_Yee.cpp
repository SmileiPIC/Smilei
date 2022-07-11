
#include "MF_Solver2D_Yee.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MF_Solver2D_Yee::MF_Solver2D_Yee( Params &params )
    : Solver2D( params )
{
    isEFilterApplied = false;
    if( params.Friedman_filter ) {
        isEFilterApplied = true;
    }
}

MF_Solver2D_Yee::~MF_Solver2D_Yee()
{
    // EMPTY
}

void MF_Solver2D_Yee::operator()( ElectroMagn *fields )
{
    const double *const __restrict__ Ex2D = isEFilterApplied ? fields->Exfilter[0]->data() :
                                                               fields->Ex_->data(); // [x * ny_p + y] : dual in x   primal in y,z
    const double *const __restrict__ Ey2D = isEFilterApplied ? fields->Eyfilter[0]->data() :
                                                               fields->Ey_->data(); // [x * ny_d + y] : dual in y   primal in x,z
    const double *const __restrict__ Ez2D = fields->Ez_->data();                    // [x * ny_p + y] : dual in z   primal in x,y
    double *const __restrict__ Bx2D       = fields->Bx_->data();                    // [x * ny_d + y] : dual in y,z primal in x
    double *const __restrict__ By2D       = fields->By_->data();                    // [x * ny_p + y] : dual in x,z primal in y
    double *const __restrict__ Bz2D       = fields->Bz_->data();                    // [x * ny_d + y] : dual in x,y primal in z

    // Magnetic field Bx^(p,d)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    const unsigned Bx_Bx2D_first = 1 - 1;
    const unsigned Bx_Bx2D_last  = ( nx_d - 1 - 1 ) * ny_d + ny_d - 1;
    const unsigned Bx_Ez2D_first = 1 - 1;
    const unsigned Bx_Ez2D_last  = ( nx_d - 1 - 1 ) * ny_p + ny_d - 1;

    #pragma omp target map( tofrom                                                \
                            : Bx2D [Bx_Bx2D_first:Bx_Bx2D_last - Bx_Bx2D_first] ) \
        map( to                                                                   \
             : Ez2D [Bx_Ez2D_first:Bx_Ez2D_last - Bx_Ez2D_first] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 0; x < nx_d - 1; ++x ) {
        for( unsigned int y = 1; y < ny_d - 1; ++y ) {
            Bx2D[x * ny_d + y] -= dt_ov_dy * ( Ez2D[x * ny_p + y] - Ez2D[x * ny_p + y - 1] );
        }
    }

    // Magnetic field By^(d,p)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    const unsigned By_By2D_first = ny_p;
    const unsigned By_By2D_last  = ( nx_d - 1 - 1 ) * ny_p + ny_p;
    const unsigned By_Ez2D_first = ny_p - ny_p;
    const unsigned By_Ez2D_last  = ( nx_d - 1 - 1 ) * ny_p + ny_p;

    #pragma omp target map( tofrom                                                \
                            : By2D [By_By2D_first:By_By2D_last - By_By2D_first] ) \
        map( to                                                                   \
             : Ez2D [By_Ez2D_first:By_Ez2D_last - By_Ez2D_first] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 1; x < nx_d - 1; ++x ) {
        for( unsigned int y = 0; y < ny_p; ++y ) {
            By2D[x * ny_p + y] += dt_ov_dx * ( Ez2D[x * ny_p + y] - Ez2D[( x - 1 ) * ny_p + y] );
        }
    }

    // Magnetic field Bz^(d,d)
#if defined( SMILEI_ACCELERATOR_GPU_OMP_PENDING )
    const unsigned Bz_Bz2D_first = ny_d + 1;
    const unsigned Bz_Bz2D_last  = ( nx_d - 1 - 1 ) * ny_d + ny_d - 1;
    const unsigned Bz_Ex2D_first = ny_p + 1 - 1;
    const unsigned Bz_Ex2D_last  = ( nx_d - 1 - 1 ) * ny_p + ny_d - 1;
    const unsigned Bz_Ey2D_first = ny_d - ny_d + 1;
    const unsigned Bz_Ey2D_last  = ( nx_d - 1 - 1 ) * ny_d + ny_d - 1;

    #pragma omp target map( tofrom                                                \
                            : Bz2D [Bz_Bz2D_first:Bz_Bz2D_last - Bz_Bz2D_first] ) \
        map( to                                                                   \
             : Ex2D [Bz_Ex2D_first:Bz_Ex2D_last - Bz_Ex2D_first],                 \
               Ey2D [Bz_Ey2D_first:Bz_Ey2D_last - Bz_Ey2D_first] )
    #pragma omp teams
    #pragma omp distribute parallel for collapse( 2 )
#endif
    for( unsigned int x = 1; x < nx_d - 1; ++x ) {
        for( unsigned int y = 1; y < ny_d - 1; ++y ) {
            Bz2D[x * ny_d + y] += dt_ov_dy * ( Ex2D[x * ny_p + y] - Ex2D[x * ny_p + y - 1] ) -
                                  dt_ov_dx * ( Ey2D[x * ny_d + y] - Ey2D[( x - 1 ) * ny_d + y] );
        }
    }
}
