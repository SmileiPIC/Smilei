
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
    cField2D* ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->El_[imode];
    cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Er_[imode];
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Et_[imode];
    cField2D* BxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bl_m[imode];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Br_m[imode];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bt_m[imode];
    cField2D* JxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jl_[imode];
    cField2D* JyRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jr_[imode];
    cField2D* JzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jt_[imode];


    // Electric field Ex^(d,p)
    for (unsigned int i=0 ; i<nx_d ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*ExRZ)(i,j) += -dt*(*JxRZ)(i,j)
                +                 dt_ov_dy * ( (*BtRZ)(i,j+1) - (*BtRZ)(i,j) );
                }
    }
    
    // Electric field Er^(p,d)
    for (unsigned int i=0 ; i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_d ; j++) {
            (*ErRZ)(i,j) += -dt*(*JyRZ)(i,j)
                -                  dt_ov_dx * ( (*BtRZ)(i+1,j) - (*BtRZ)(i,j) );
        }
    }
    
    // Electric field Et^(p,p)
    for (unsigned int i=0 ;  i<nx_p ; i++) {
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*EtRZ)(i,j) += -dt*(*JzRZ)(i,j)
                +                  dt_ov_dx * ( (*BrRZ)(i+1,j) - (*BrRZ)(i,j) )
                -                  dt_ov_dy * ( (*BxRZ)(i,j+1) - (*BxRZ)(i,j) );
        }
    }

}

