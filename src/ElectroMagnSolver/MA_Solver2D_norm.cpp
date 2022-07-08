
#include "MA_Solver2D_norm.h"

#include "ElectroMagn.h"
#include "Field2D.h"

MA_Solver2D_norm::MA_Solver2D_norm( Params &params )
    : Solver2D( params )
{
}

MA_Solver2D_norm::~MA_Solver2D_norm()
{
}

void MA_Solver2D_norm::operator()( ElectroMagn *fields )
{
    double *const __restrict__ Ex2D       = fields->Ex_->data(); // [x * ny_p + y] : dual in x   primal in y,z
    double *const __restrict__ Ey2D       = fields->Ey_->data(); // [x * ny_d + y] : dual in y   primal in x,z
    double *const __restrict__ Ez2D       = fields->Ez_->data(); // [x * nz_p + y] : dual in z   primal in x,y
    const double *const __restrict__ Bx2D = fields->Bx_->data(); // [x * ny_d + y] : dual in y,z primal in x
    const double *const __restrict__ By2D = fields->By_->data(); // [x * ny_p + y] : dual in x,z primal in y
    const double *const __restrict__ Bz2D = fields->Bz_->data(); // [x * ny_d + y] : dual in x,y primal in z
    const double *const __restrict__ Jx2D = fields->Jx_->data(); // [x * ny_p + y] : dual in x   primal in y,z
    const double *const __restrict__ Jy2D = fields->Jy_->data(); // [x * ny_d + y] : dual in y   primal in x,z
    const double *const __restrict__ Jz2D = fields->Jz_->data(); // [x * ny_p + y] : dual in z   primal in x,y

    // Electric field Ex^(d,p)
    for( unsigned int x = 0; x < nx_d; ++x ) {
        for( unsigned int y = 0; y < ny_p; ++y ) {
            Ex2D[x * ny_p + y] += -dt * Jx2D[x * ny_p + y] + dt_ov_dy * ( Bz2D[x * ny_d + y + 1] - Bz2D[x * ny_d + y] );
        }
    }

    // Electric field Ey^(p,d)
    for( unsigned int x = 0; x < nx_p; ++x ) {
        for( unsigned int y = 0; y < ny_d; ++y ) {
            Ey2D[x * ny_d + y] += -dt * Jy2D[x * ny_d + y] - dt_ov_dx * ( Bz2D[( x + 1 ) * ny_d + y] - Bz2D[x * ny_d + y] );
        }
    }

    // Electric field Ez^(p,p)
    for( unsigned int x = 0; x < nx_p; ++x ) {
        for( unsigned int y = 0; y < ny_p; ++y ) {
            Ez2D[x * ny_p + y] += -dt * Jz2D[x * ny_p + y] +
                                  dt_ov_dx * ( By2D[( x + 1 ) * ny_p + y] - By2D[x * ny_p + y] ) -
                                  dt_ov_dy * ( Bx2D[x * ny_d + y + 1] - Bx2D[x * ny_d + y] );
        }
    }
}
