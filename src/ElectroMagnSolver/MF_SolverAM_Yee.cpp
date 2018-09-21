
#include "MF_SolverAM_Yee.h"

#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"

MF_SolverAM_Yee::MF_SolverAM_Yee(Params &params)
: SolverAM(params)
{
    isEFilterApplied = false;
    if (params.Friedman_filter)
        ERROR("Filters are not available yet");;
}

MF_SolverAM_Yee::~MF_SolverAM_Yee()
{
}

void MF_SolverAM_Yee::operator() ( ElectroMagn* fields )
{

    for (unsigned int imode=0 ; imode<Nmode ; imode++) {

        // Static-cast of the fields

        cField2D* ElAM = (static_cast<ElectroMagn3DAM*>(fields))->El_[imode];
        cField2D* ErAM = (static_cast<ElectroMagn3DAM*>(fields))->Er_[imode];
        cField2D* EtAM = (static_cast<ElectroMagn3DAM*>(fields))->Et_[imode];
        cField2D* BlAM = (static_cast<ElectroMagn3DAM*>(fields))->Bl_[imode];
        cField2D* BrAM = (static_cast<ElectroMagn3DAM*>(fields))->Br_[imode];
        cField2D* BtAM = (static_cast<ElectroMagn3DAM*>(fields))->Bt_[imode];
        int  j_glob = (static_cast<ElectroMagn3DAM*>(fields))->j_glob_;
        bool isYmin = (static_cast<ElectroMagn3DAM*>(fields))->isYmin;
        bool isXmin = (static_cast<ElectroMagn3DAM*>(fields))->isXmin;
        bool isYmax = (static_cast<ElectroMagn3DAM*>(fields))->isYmax;
        bool isXmax = (static_cast<ElectroMagn3DAM*>(fields))->isXmax;

        // Magnetic field Bx^(p,d)
        for (unsigned int i=0 ; i<nl_p;  i++) {
            #pragma omp simd
            for (unsigned int j=1+isYmin*2 ; j<nr_d-1 ; j++) {
                (*BlAM)(i,j) += - dt/((j_glob+j-0.5)*dr) * ( (double)(j+j_glob)*(*EtAM)(i,j) - (double)(j+j_glob-1.)*(*EtAM)(i,j-1) + Icpx*(double)imode*(*ErAM)(i,j));
            }
        }
            
            // Magnetic field Br^(d,p)
        for (unsigned int i=1 ; i<nl_d-1 ; i++) {
            #pragma omp simd
            for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {  //Specific condition on axis
                (*BrAM)(i,j) += dt_ov_dl * ( (*EtAM)(i,j) - (*EtAM)(i-1,j) )
                 +Icpx*dt*(double)imode/((double)(j_glob+j)*dr)*(*ElAM)(i,j) ;
            }
        }
            // Magnetic field Bt^(d,d)
        for (unsigned int i=1 ; i<nl_d-1 ; i++) {
            #pragma omp simd
            for (unsigned int j=1 + isYmin*2 ; j<nr_d-1 ; j++) {
                (*BtAM)(i,j) += dt_ov_dr * ( (*ElAM)(i,j) - (*ElAM)(i,j-1) )
                -dt_ov_dl * ( (*ErAM)(i,j) - (*ErAM)(i-1,j) );
            }
        }

        // On axis conditions
        if (isYmin){
            unsigned int j=2;
            if (imode==0){
            	for (unsigned int i=1 ; i<nl_d-1 ; i++) {
            		(*BrAM)(i,j)=0;
            	}
            	for (unsigned int i=1 ; i<nl_d-1 ; i++) {
            		(*BtAM)(i,j)= -(*BtAM)(i,j+1);
            	}
            	for (unsigned int i=0 ; i<nl_p ; i++) {
            		(*BlAM)(i,j)= (*BlAM)(i,j+1);
            	}
            }
            
            else if (imode==1){
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*BlAM)(i,j)= -(*BlAM)(i,j+1);
            	}
            
            	for (unsigned int i=1 ; i<nl_d-1 ; i++) {
            		(*BrAM)(i,j)+=  Icpx*dt_ov_dr*(*ElAM)(i,j+1)
            		+			dt_ov_dl*((*EtAM)(i,j)-(*EtAM)(i-1,j));
            	}
            	for (unsigned int i=1; i<nl_d-1 ; i++) {
            		(*BtAM)(i,j)= -2.*Icpx*(*BrAM)(i,j)-(*BtAM)(i,j+1);
            	}	
            
            }
            else {  // modes > 1
            	for (unsigned int  i=0 ; i<nl_p; i++) {
            		(*BlAM)(i,j)= -(*BlAM)(i,j+1);
            	}
            	for (unsigned int i=1 ; i<nl_d-1; i++) {
            		(*BrAM)(i,j)= 0;
            	}
            	for (unsigned int  i=1 ; i<nl_d-1 ; i++) {
            		(*BtAM)(i,j)= - (*BtAM)(i,j+1);
            	}	
            }
        }
    }
}

