
#include "MF_SolverRZ_Yee.h"

#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverRZ_Yee::MF_SolverRZ_Yee(Params &params)
: SolverRZ(params)
{
    isEFilterApplied = false;
    if (params.Friedman_filter)
        ERROR("Filters are not available yet");;
}

MF_SolverRZ_Yee::~MF_SolverRZ_Yee()
{
}

void MF_SolverRZ_Yee::operator() ( ElectroMagn* fields )
{

    #ifdef _TODO_RZ
    #endif
    // Loop on modes ?
    for (unsigned int imode=0 ; imode<Nmode ; imode++) {
    

    //int imode = 0;

    // Static-cast of the fields

    cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(fields))->El_[imode];
    cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Er_[imode];
    cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Et_[imode];
    cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bl_[imode];
    cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Br_[imode];
    cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bt_[imode];
    int     j_glob = (static_cast<ElectroMagn3DRZ*>(fields))->j_glob_;
	bool isYmin = (static_cast<ElectroMagn3DRZ*>(fields))->isYmin;

    //    #pragma omp simd
    //    for (unsigned int j=0 ; j<ny_d-1 ; j++) {
    //        (*BlRZ)(nl_d,j) -= dt_ov_dy * ( (*EtRZ)(nl_d,j) - (*EtRZ)(nl_d,j-1) );
    //    }
    // Magnetic field Bx^(p,d)
    for (unsigned int i=0 ; i<nl_p;  i++) {
        #pragma omp simd
        for (unsigned int j=isYmin*3 ; j<nr_d ; j++) {
            (*BlRZ)(i,j) += - dt/((j_glob+j-0.5)*dr) * ( (double)(j)*(*EtRZ)(i,j) - (double)(j-1)*(*EtRZ)(i,j-1) )
             -Icpx*dt*(double)imode/((j_glob+j-0.5)*dr)*(*ErRZ)(i,j);
			 if (std::abs((*BlRZ)(i,j))>1.){
                MESSAGE("BlRZMF");                
                MESSAGE(i);
                MESSAGE(j);    
                MESSAGE((*BlRZ)(i,j));
                }

        }
    }
        
        // Magnetic field Br^(d,p)
    for (unsigned int i=1 ; i<nl_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {
            (*BrRZ)(i,j) += dt_ov_dl * ( (*EtRZ)(i,j) - (*EtRZ)(i-1,j) )
             +Icpx*dt*(double)imode/((double)(j_glob+j)*dr)*(*ElRZ)(i,j) ;
			 if (std::abs((*BrRZ)(i,j))>1.){
                MESSAGE("BrRZMF");                
                MESSAGE(i);
                MESSAGE(j);    
                MESSAGE((*BrRZ)(i,j));
                }

        }
    }
        // Magnetic field Bt^(d,d)
    for (unsigned int i=1 ; i<nl_d-1 ; i++) {
        #pragma omp simd
        for (unsigned int j=isYmin*3 ; j<nr_d-1 ; j++) {
            (*BtRZ)(i,j) += dt_ov_dr * ( (*ElRZ)(i,j) - (*ElRZ)(i,j-1) )
            -dt_ov_dl * ( (*ErRZ)(i,j) - (*ErRZ)(i-1,j) );
		    if (std::abs((*BtRZ)(i,j))>1.){
                MESSAGE("BtRZMF");                
                MESSAGE(i);
                MESSAGE(j);    
                MESSAGE((*BtRZ)(i,j));
                }
        }
    }

    
    }// end parallel
}

