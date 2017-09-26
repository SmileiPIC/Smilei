
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
    cField2D* ErRZ;
    //if (!isEFilterApplied) {
        ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Ex_[imode];
        ErRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Er_[imode];
    //} else {
    //    ExRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Exfilter[0];
    //    ErRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Erfilter[0];
    //}
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Et_[imode];
    cField2D* BxRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bx_[imode];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Br_[imode];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bt_[imode];
    
    // Magnetic field Bx^(p,d)
    //cout << "nx_p,nx_d-1" << nx_p << " " << nx_d-1 ;
    //{
    {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BxRZ)(0,j) -= dt_ov_dy * ( (*EtRZ)(0,j) - (*EtRZ)(0,j-1) );
        }
    }
    //    for (unsigned int i=0 ; i<nx_p;  i++) {
    for (unsigned int i=1 ; i<nx_d-1;  i++) {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BxRZ)(i,j) -= dt_ov_dy * ( (*EtRZ)(i,j) - (*EtRZ)(i,j-1) );
        }
        //    }
        
        // Magnetic field Br^(d,p)
        //    for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=0 ; j<ny_p ; j++) {
            (*BrRZ)(i,j) += dt_ov_dx * ( (*EtRZ)(i,j) - (*EtRZ)(i-1,j) );
        }
        //}
        
        // Magnetic field Bt^(d,d)
        //for (unsigned int i=1 ; i<nx_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=1 ; j<ny_d-1 ; j++) {
            (*BtRZ)(i,j) += dt_ov_dy * ( (*ExRZ)(i,j) - (*ExRZ)(i,j-1) )
            -               dt_ov_dx * ( (*ErRZ)(i,j) - (*ErRZ)(i-1,j) );
        }
    }
    //}// end parallel
}

