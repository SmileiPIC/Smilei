
#include "MA_SolverRZ_norm.h"

#include "ElectroMagn3DRZ.h"
#include "cField2D.h"

MA_SolverRZ_norm::MA_SolverRZ_norm(Params &params)
: SolverRZ(params)
{
}

MA_SolverRZ_norm::~MA_SolverRZ_norm()
{
}

void MA_SolverRZ_norm::operator() ( ElectroMagn* fields )
{

    #ifdef _TODO_RZ
    // Loop on modes ?
    #endif

    int imode = 0;

    // Static-cast of the fields
    cField2D* ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ex_RZ_[imode];
    cField2D* EyRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ey_RZ_[imode];
    cField2D* EzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ez_RZ_[imode];
    cField2D* BxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bx_RZ_m[imode];
    cField2D* ByRZ = (static_cast<ElectroMagn3DRZ*>(fields))->By_RZ_m[imode];
    cField2D* BzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bz_RZ_m[imode];
    cField2D* JxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jx_RZ_[imode];
    cField2D* JyRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jy_RZ_[imode];
    cField2D* JzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jz_RZ_[imode];


    // Electric field Ex^(d,p,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*ExRZ)(i,j) += -dt*(*JxRZ)(i,j)
                +                 dt_ov_dy * ( (*BzRZ)(i,j+1) - (*BzRZ)(i,j) );
                }
    }
    
    // Electric field Ey^(p,d,p)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*EyRZ)(i,j) += -dt*(*JyRZ)(i,j)
                -                  dt_ov_dx * ( (*BzRZ)(i+1,j) - (*BzRZ)(i,j) );
        }
    }
    
    // Electric field Ez^(p,p,d)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*EzRZ)(i,j) += -dt*(*JzRZ)(i,j)
                +                  dt_ov_dx * ( (*ByRZ)(i+1,j) - (*ByRZ)(i,j) )
                -                  dt_ov_dy * ( (*BxRZ)(i,j+1) - (*BxRZ)(i,j) );
        }
    }

}

