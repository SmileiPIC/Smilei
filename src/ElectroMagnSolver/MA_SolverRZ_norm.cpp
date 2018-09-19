
#include "MA_SolverRZ_norm.h"
#include "ElectroMagn3DRZ.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"
#include "Patch.h"
MA_SolverRZ_norm::MA_SolverRZ_norm(Params &params)
: SolverRZ(params)
{
}

MA_SolverRZ_norm::~MA_SolverRZ_norm()
{
}

void MA_SolverRZ_norm::operator() ( ElectroMagn* fields )
{
    for (unsigned int imode=0 ; imode<Nmode ; imode++) {
     
        // Static-cast of the fields_SolverRZ_norm.cpp
        cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(fields))->El_[imode];
        cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Er_[imode];
        cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Et_[imode];
        cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bl_[imode];
        cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Br_[imode];
        cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Bt_[imode];
        cField2D* JlRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jl_[imode];
        cField2D* JrRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jr_[imode];
        cField2D* JtRZ = (static_cast<ElectroMagn3DRZ*>(fields))->Jt_[imode];
        int j_glob    = (static_cast<ElectroMagn3DRZ*>(fields))->j_glob_;
        bool isYmin = (static_cast<ElectroMagn3DRZ*>(fields))->isYmin;
        bool isXmin = (static_cast<ElectroMagn3DRZ*>(fields))->isXmin;
        bool isXmax = (static_cast<ElectroMagn3DRZ*>(fields))->isXmax;
        bool isYmax = (static_cast<ElectroMagn3DRZ*>(fields))->isYmax;

        // Electric field Elr^(d,p)
        for (unsigned int i=0 ; i<nl_d ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {  
                (*ElRZ)(i,j) += -dt*(*JlRZ)(i,j) 
                    +                 dt/((j_glob+j)*dr)*((j+j_glob+0.5)*(*BtRZ)(i,j+1) - (j+j_glob-0.5)*(*BtRZ)(i,j) )
                    +                 Icpx*dt*(double)imode/((j_glob+j)*dr)*(*BrRZ)(i,j);
            }
        }
        for (unsigned int i=0 ; i<nl_p ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_d ; j++) {
                (*ErRZ)(i,j) += -dt*(*JrRZ)(i,j)
                    -                  dt_ov_dl * ( (*BtRZ)(i+1,j) - (*BtRZ)(i,j) )
                    -                  Icpx*dt*(double)imode/((j_glob+j-0.5)*dr)* (*BlRZ)(i,j);

            }
        }
        for (unsigned int i=0 ;  i<nl_p ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {
                (*EtRZ)(i,j) += -dt*(*JtRZ)(i,j)
                    +                  dt_ov_dl * ( (*BrRZ)(i+1,j) - (*BrRZ)(i,j) )
                    -                  dt_ov_dr * ( (*BlRZ)(i,j+1) - (*BlRZ)(i,j) );        
            }
        }
        if (isYmin){ // Conditions on axis
            unsigned int j=2;
            if (imode==0){
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*EtRZ)(i,j)=0;
            	}
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*ErRZ)(i,j)= -(*ErRZ)(i,j+1);
            	}
            	for (unsigned int i=0 ; i<nl_d ; i++) {
            		(*ElRZ)(i,j)+= 4.*dt_ov_dr*(*BtRZ)(i,j+1)-dt*(*JlRZ)(i,j);
            	}
            }
            else if (imode==1){
            	for (unsigned int i=0 ; i<nl_d  ; i++) {
            		(*ElRZ)(i,j)= 0;
            	}
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*EtRZ)(i,j)= -1./3.*(4.*Icpx*(*ErRZ)(i,j+1)+(*EtRZ)(i,j+1));
            	}
            	for (unsigned int i=0 ; i<nl_p ; i++) {
            		(*ErRZ)(i,j)=2.*Icpx*(*EtRZ)(i,j)-(*ErRZ)(i,j+1);
            	}
            }
            else {	
            	for (unsigned int  i=0 ; i<nl_d; i++) {
            		(*ElRZ)(i,j)= 0;
            	}
            	for (unsigned int  i=0 ; i<nl_p; i++) {
            		(*ErRZ)(i,j)= -(*ErRZ)(i,j+1);
            	}
            	for (unsigned int i=0 ; i<nl_p; i++) {
            		(*EtRZ)(i,j)= 0;
            	}
            }
        } 
    }
}

