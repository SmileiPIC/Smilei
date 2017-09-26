
#include "MF_SolverRZ_Yee.h"

#include "ElectroMagn3DRZ.h"
#include "cField2D.h"

MF_SolverRZ_Yee::MF_SolverRZ_Yee(Params &params)
: SolverRZ(params)
{
    isEFilterApplied = false;
    if (params.Friedman_filter)
        isEFilterApplied = true;
}

MF_SolverRZ_Yee::~MF_SolverRZ_Yee()
{
}

void MF_SolverRZ_Yee::operator() ( ElectroMagn* fields )
{

    #ifdef _TODO_RZ
    // Loop on modes ?
    #endif

    int imode = 0;

    // Static-cast of the fields
    cField2D* ExRZ;
    cField2D* EyRZ;
    //if (!isEFilterApplied) {
        ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ex_RZ_[imode];
        EyRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ey_RZ_[imode];
    //} else {
    //    ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Exfilter[0];
    //    EyRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Eyfilter[0];
    //}
    cField2D* EzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ez_RZ_[imode];
    cField2D* BxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bx_RZ_[imode];
    cField2D* ByRZ = (static_cast<ElectroMagn3DRZ*>(fields))->By_RZ_[imode];
    cField2D* BzRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bz_RZ_[imode];
    
    // Magnetic field Bx^(p,d)
    //cout << "nx_p,nx_d-1" << nx_p << " " << nx_d-1 ;
    //{
    {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BxRZ)(0,j) -= dt_ov_dy * ( (*EzRZ)(0,j) - (*EzRZ)(0,j-1) );
        }
    }
    //    for (unsigned int i=0 ; i<nx_p;  i++) {
    for (unsigned int i=1 ; i<nx_d-1;  i++) {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BxRZ)(i,j) -= dt_ov_dy * ( (*EzRZ)(i,j) - (*EzRZ)(i,j-1) );
        }
        //    }
        
        // Magnetic field By^(d,p)
        //    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*ByRZ)(i,j) += dt_ov_dx * ( (*EzRZ)(i,j) - (*EzRZ)(i-1,j) );
        }
        //}
        
        // Magnetic field Bz^(d,d)
        //for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BzRZ)(i,j) += dt_ov_dy * ( (*ExRZ)(i,j) - (*ExRZ)(i,j-1) )
            -               dt_ov_dx * ( (*EyRZ)(i,j) - (*EyRZ)(i-1,j) );
        }
    }
    //}// end parallel
}

