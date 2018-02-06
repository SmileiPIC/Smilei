#include "ElectroMagnBCRZ_BM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"
#include "Tools.h"
#include "Laser.h"
#include <complex>
#include "dcomplex.h"
using namespace std;

ElectroMagnBCRZ_BM::ElectroMagnBCRZ_BM( Params &params, Patch* patch, unsigned int _min_max )
: ElectroMagnBC( params, patch, _min_max )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
    // number of nodes of the primal and dual grid in the x-direction
    nl_p = params.n_space[0]+1+2*params.oversize[0];
    nl_d = nl_p+1;
    // number of nodes of the primal and dual grid in the y-direction
    nr_p = params.n_space[1]+1+2*params.oversize[1];
    nr_d = nr_p+1;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dl       = params.cell_length[0];
    dt_ov_dl = dt/dl;
    dl_ov_dt = 1.0/dt_ov_dl;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dr       = params.cell_length[1];
    dt_ov_dr = dt/dr;
    dr_ov_dt = 1.0/dt_ov_dr;
    //Number of modes
	Nmode=params.Nmode;
    if (min_max == 0 && patch->isXmin() ) {
        // BCs at the x-border min
        Bl_val.resize(nr_d,0.); // dual in the y-direction
        Br_val.resize(nr_p,0.); // primal in the y-direction
        Bt_val.resize(nr_d,0.); // dual in the y-direction
    }
    else if (min_max == 1 && patch->isXmax() ) {
        // BCs at the x-border max
        Bl_val.resize(nr_d,0.);
        Br_val.resize(nr_p,0.);
        Bt_val.resize(nr_d,0.);
    }
    else if (min_max == 2 && patch->isYmin() ) {
        // BCs in the y-border min
        Bl_val.resize(nl_p,0.); // primal in the x-direction
        Br_val.resize(nl_d,0.); // dual in the x-direction
        Bt_val.resize(nl_d,0.); // dual in the x-direction
    }
    else if (min_max == 3 && patch->isYmax() ) {
        // BCs in the y-border max
        Bl_val.resize(nl_p,0.);
        Br_val.resize(nl_d,0.);
        Bt_val.resize(nl_d,0.);
    }
    
    
    
    
    // -----------------------------------------------------
    // Parameters for the Silver-Mueller boundary conditions
    // -----------------------------------------------------

    #ifdef _TODO_RZ
    #endif
    // Rmax boundary
    double theta  = 0.0;
	double phi    =  0.0;
	CB_BM  = cos(theta)*cos(phi)/(1.0 + cos(theta)*cos(phi));
	CE_BM  = 1.0 - CB_BM;
    Alpha_BM_Rmax    = (1. - dt_ov_dr)/(1. + dt_ov_dr);
    Beta_BM_Rmax     = dt/(1.+ dt_ov_dr);
    Gamma_BM_Rmax    =  dt_ov_dl/(1. + dt_ov_dr);
    
}



void ElectroMagnBCRZ_BM::save_fields(Field* my_field, Patch* patch) {
    cField2D* field2D=static_cast<cField2D*>(my_field);
    
    if (min_max == 0 && patch->isXmin() ) {
        if (field2D->name=="Bl"){
            for (unsigned int j=0; j<nr_d; j++) {
                Bl_val[j]=(*field2D)(0,j);
            }
        }
        
        if (field2D->name=="Br"){
            for (unsigned int j=0; j<nr_p; j++) {
                Br_val[j]=(*field2D)(0,j);
            }
        }
        
        if (field2D->name=="Bt"){
            for (unsigned int j=0; j<nr_d; j++) {
                Bt_val[j]=(*field2D)(0,j);
            }
        }
    }
    else if (min_max == 1 && patch->isXmax() ) {
        if (field2D->name=="Bl"){
            for (unsigned int j=0; j<nr_d; j++) {
                Bl_val[j]=(*field2D)(nl_p-1,j);
            }
        }
        
        if (field2D->name=="Br"){
            for (unsigned int j=0; j<nr_p; j++) {
                Br_val[j]=(*field2D)(nl_d-1,j);
            }
        }
        
        if (field2D->name=="Bt"){
            for (unsigned int j=0; j<nr_d; j++) {
                Bt_val[j]=(*field2D)(nl_d-1,j);
            }
        }
    }
    else if (min_max == 2 && patch->isYmin() ) {
        if (field2D->name=="Bl"){
            for (unsigned int i=0; i<nl_p; i++) {
                Bl_val[i]=(*field2D)(i,0);
            }
        }
        
        if (field2D->name=="Br"){
            for (unsigned int i=0; i<nl_d; i++) {
                Br_val[i]=(*field2D)(i,0);
            }
        }
        
        if (field2D->name=="Bt"){
            for (unsigned int i=0; i<nl_d; i++) {
                Bt_val[i]=(*field2D)(i,0);
            }
        }
    }
    else if (min_max == 3 && patch->isYmax() ) {
        if (field2D->name=="Bl"){
            for (unsigned int i=0; i<nl_p; i++) {
                Bl_val[i]=(*field2D)(i,nr_d-1);
            }
        }
        
        if (field2D->name=="Br"){
            for (unsigned int i=0; i<nl_d; i++) {
                Br_val[i]=(*field2D)(i,nr_p-1);
            }
        }
        
        if (field2D->name=="Bt"){
            for (unsigned int i=0; i<nl_d; i++) {
                Bt_val[i]=(*field2D)(i,nr_d-1);
            }
        }
    }
}

void ElectroMagnBCRZ_BM::disableExternalFields()
{
    Bl_val.resize(0);
    Br_val.resize(0);
    Bt_val.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCRZ_BM::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    // Loop on imode 
    for (unsigned int imode=0 ; imode<Nmode ; imode++) {

		// Static cast of the fields
                //cField2D* ElRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[imode];
		cField2D* ErRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[imode];
		cField2D* EtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[imode];
		cField2D* BlRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_[imode];
		cField2D* BrRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_[imode];
		cField2D* BtRZ = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_[imode];
		cField2D* BtRZ_old = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_m[imode]; 
		cField2D* BlRZ_old = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_m[imode];
		cField2D* BrRZ_old = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_m[imode];
		int     j_glob = (static_cast<ElectroMagn3DRZ*>(EMfields))->j_glob_;
		bool isXmin = (static_cast<ElectroMagn3DRZ*>(EMfields))->isXmin;
		bool isXmax = (static_cast<ElectroMagn3DRZ*>(EMfields))->isXmax;
		if (min_max == 3 && patch->isYmax() ) {
		    
		    // for Bl^(p,d)
			unsigned int j= nr_d-2;
		    for (unsigned int i=isXmin ; i<nl_p-isXmax; i++) {
		        /*(*Bl2D)(i,nr_d-1) = -Alpha_SM_N   * (*Et2D)(i,nr_p-1)
		         +                    Beta_SM_N    * (*Bl2D)(i,nr_d-2)
		         +                    Delta_SM_N   * (*Br2D)(i+1,nr_p-1)
		         +                    Epsilon_SM_N * (*Br2D)(i,nr_p-1);*/
		        (*BlRZ)(i,j+1) = (*BlRZ_old)(i,j) - Alpha_BM_Rmax   * ((*BlRZ)(i,j)-(*BlRZ_old)(i,j+1))
		        -                   Gamma_BM_Rmax*CB_BM    *( (*BrRZ)(i+1,j) + (*BrRZ_old)(i+1,j) - (*BrRZ)(i,j) - (*BrRZ_old)(i,j))
		        -                   Beta_BM_Rmax*Icpx*(double)imode*CE_BM/((j_glob+j+0.5)*dr) *( (*ErRZ)(i,j+1)+(*ErRZ)(i,j))
		        -                   2.*CE_BM*dt/((j_glob+j+0.5)*dr)*(*EtRZ)(i,j);
                if (std::abs((*BlRZ)(i,j))>1.)
                {
                	MESSAGE("BlRZBM");                
                	MESSAGE(i);    
                	MESSAGE(j);
                	MESSAGE((*BlRZ)(i,j));
                }
		    }//i  ---end Bl
		    
		    // for Bt^(d,d)
			
		    for (unsigned int i=1 ; i<nl_d-2 ; i++) {
		        /*(*Bt2D)(i,nr_d-1) = Alpha_SM_N * (*El2D)(i,nr_p-1)
		         +                   Beta_SM_N  * (*Bt2D)(i,nr_d-2);*/
		        (*BtRZ)(i,j+1) = (*BtRZ_old)(i,j)- Alpha_BM_Rmax * ((*BtRZ)(i,j) - (*BtRZ_old)(i,j+1))
		        -                   Icpx*(double)imode*CB_BM*Beta_BM_Rmax/((j_glob+j+0.5)*dr)  *((*BrRZ)(i,j) - (*BrRZ_old)(i,j) )
		        -                   CE_BM*Gamma_BM_Rmax*((*ErRZ)(i,j+1)+(*ErRZ)(i,j)-(*ErRZ)(i-1,j) -(*ErRZ)(i-1,j+1) )
				-                   CB_BM* Beta_BM_Rmax/((j_glob+j+0.5)*dr +(j_glob+j-0.5)*dr)*((*BtRZ)(i,j+1) + (*BtRZ_old)(i,j+1) 				+					(*BtRZ)(i,j) + (*BtRZ_old)(i,j)) ;
                if (std::abs((*BtRZ)(i,j))>1.)
                {
                	MESSAGE("BtRZBM");                
                	MESSAGE(i);    
                	MESSAGE(j);
                	MESSAGE((*BtRZ)(i,j));
                }
		    }//i  ---end Bt
		    
		}
		else  {
		    ERROR( "No Buneman along the axis" );
		}
	}	
}
