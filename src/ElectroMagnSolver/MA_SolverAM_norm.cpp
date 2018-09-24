
#include "MA_SolverAM_norm.h"
#include "ElectroMagnAM.h"
#include "cField2D.h"
#include <complex>
#include "dcomplex.h"
#include "Patch.h"
MA_SolverAM_norm::MA_SolverAM_norm(Params &params)
: SolverAM(params)
{
}

MA_SolverAM_norm::~MA_SolverAM_norm()
{
}

void MA_SolverAM_norm::operator() ( ElectroMagn* fields )
{
    for (unsigned int imode=0 ; imode<Nmode ; imode++) {
     
        // Static-cast of the fields_SolverAM_norm.cpp
        cField2D* ElAM = (static_cast<ElectroMagnAM*>(fields))->El_[imode];
        cField2D* ErAM = (static_cast<ElectroMagnAM*>(fields))->Er_[imode];
        cField2D* EtAM = (static_cast<ElectroMagnAM*>(fields))->Et_[imode];
        cField2D* BlAM = (static_cast<ElectroMagnAM*>(fields))->Bl_[imode];
        cField2D* BrAM = (static_cast<ElectroMagnAM*>(fields))->Br_[imode];
        cField2D* BtAM = (static_cast<ElectroMagnAM*>(fields))->Bt_[imode];
        cField2D* JlAM = (static_cast<ElectroMagnAM*>(fields))->Jl_[imode];
        cField2D* JrAM = (static_cast<ElectroMagnAM*>(fields))->Jr_[imode];
        cField2D* JtAM = (static_cast<ElectroMagnAM*>(fields))->Jt_[imode];
        int j_glob    = (static_cast<ElectroMagnAM*>(fields))->j_glob_;
        bool isYmin = (static_cast<ElectroMagnAM*>(fields))->isYmin;
        //bool isXmin = (static_cast<ElectroMagnAM*>(fields))->isXmin;
        //bool isXmax = (static_cast<ElectroMagnAM*>(fields))->isXmax;
        //bool isYmax = (static_cast<ElectroMagnAM*>(fields))->isYmax;

        // Electric field Elr^(d,p)
        for (unsigned int i=0 ; i<nl_d ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {  
                (*ElAM)(i,j) += -dt*(*JlAM)(i,j) 
                    +                 dt/((j_glob+j)*dr)*((j+j_glob+0.5)*(*BtAM)(i,j+1) - (j+j_glob-0.5)*(*BtAM)(i,j) )
                    +                 Icpx*dt*(double)imode/((j_glob+j)*dr)*(*BrAM)(i,j);
            }
        }
        for (unsigned int i=0 ; i<nl_p ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_d ; j++) {
                (*ErAM)(i,j) += -dt*(*JrAM)(i,j)
                    -                  dt_ov_dl * ( (*BtAM)(i+1,j) - (*BtAM)(i,j) )
                    -                  Icpx*dt*(double)imode/((j_glob+j-0.5)*dr)* (*BlAM)(i,j);

            }
        }
        for (unsigned int i=0 ;  i<nl_p ; i++) {
            for (unsigned int j=isYmin*3 ; j<nr_p ; j++) {
                (*EtAM)(i,j) += -dt*(*JtAM)(i,j)
                    +                  dt_ov_dl * ( (*BrAM)(i+1,j) - (*BrAM)(i,j) )
                    -                  dt_ov_dr * ( (*BlAM)(i,j+1) - (*BlAM)(i,j) );        
            }
        }
        if (isYmin){ // Conditions on axis
            unsigned int j=2;
            if (imode==0){
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*EtAM)(i,j)=0;
            	}
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*ErAM)(i,j)= -(*ErAM)(i,j+1);
            	}
            	for (unsigned int i=0 ; i<nl_d ; i++) {
            		(*ElAM)(i,j)+= 4.*dt_ov_dr*(*BtAM)(i,j+1)-dt*(*JlAM)(i,j);
            	}
            }
            else if (imode==1){
            	for (unsigned int i=0 ; i<nl_d  ; i++) {
            		(*ElAM)(i,j)= 0;
            	}
            	for (unsigned int i=0 ; i<nl_p  ; i++) {
            		(*EtAM)(i,j)= -1./3.*(4.*Icpx*(*ErAM)(i,j+1)+(*EtAM)(i,j+1));
            	}
            	for (unsigned int i=0 ; i<nl_p ; i++) {
            		(*ErAM)(i,j)=2.*Icpx*(*EtAM)(i,j)-(*ErAM)(i,j+1);
            	}
            }
            else {	
            	for (unsigned int  i=0 ; i<nl_d; i++) {
            		(*ElAM)(i,j)= 0;
            	}
            	for (unsigned int  i=0 ; i<nl_p; i++) {
            		(*ErAM)(i,j)= -(*ErAM)(i,j+1);
            	}
            	for (unsigned int i=0 ; i<nl_p; i++) {
            		(*EtAM)(i,j)= 0;
            	}
            }
        } 
    }
}

