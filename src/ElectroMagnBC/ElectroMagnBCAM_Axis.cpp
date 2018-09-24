#include "ElectroMagnBCAM_Axis.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"
#include "Tools.h"
#include <complex>
#include "dcomplex.h"

using namespace std;

ElectroMagnBCAM_Axis::ElectroMagnBCAM_Axis( Params &params, Patch* patch, unsigned int _min_max )
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
	Nmode= params.nmodes;
    
   // if (min_max == 2 && patch->isYmin() ) {
        // BCs in the y-border min
   //     Bl_val.resize(nl_p,0.); // primal in the x-direction
   //     Br_val.resize(nl_d,0.); // dual in the x-direction
   //     Bt_val.resize(nl_d,0.); // dual in the x-direction
   //}
    
    #ifdef _TODO_AM
    #endif
}


void ElectroMagnBCAM_Axis::save_fields(Field* my_field, Patch* patch)
{
    ERROR( "Impossible" );
}

void ElectroMagnBCAM_Axis::disableExternalFields()
{
    ERROR( "Impossible" );
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_Axis::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{	return;
    // Loop on imode 
    for (unsigned int imode=0 ; imode< Nmode ; imode++){

    // Static cast of the fields
    cField2D* ElAM = (static_cast<ElectroMagnAM*>(EMfields))->El_[imode];
    //cField2D* ErAM = (static_cast<ElectroMagnAM*>(EMfields))->Er_[imode];
    cField2D* EtAM = (static_cast<ElectroMagnAM*>(EMfields))->Et_[imode];
    cField2D* BlAM = (static_cast<ElectroMagnAM*>(EMfields))->Bl_[imode];
    cField2D* BrAM = (static_cast<ElectroMagnAM*>(EMfields))->Br_[imode];
    cField2D* BtAM = (static_cast<ElectroMagnAM*>(EMfields))->Bt_[imode];
    //cField2D* BlAM_old = (static_cast<ElectroMagnAM*>(EMfields))->Bl_m[imode];
    //cField2D* BtAM_old = (static_cast<ElectroMagnAM*>(EMfields))->Bt_m[imode]; 
    //cField2D* JlAM = (static_cast<ElectroMagnAM*>(EMfields))->Jl_[imode];
    //cField2D* JrAM = (static_cast<ElectroMagnAM*>(EMfields))->Jr_[imode];
    //cField2D* JtAM = (static_cast<ElectroMagnAM*>(EMfields))->Jt_[imode];
    bool isXmin = (static_cast<ElectroMagnAM*>(EMfields))->isXmin;
    //bool isXmax = (static_cast<ElectroMagnAM*>(EMfields))->isXmax;

	if (min_max == 2 && patch->isYmin()){
		unsigned int j=2;
		if (imode==0){
			//MF_Solver_Yee
			for (unsigned int i=isXmin ; i<nl_d ; i++) {
				(*BrAM)(i,j)=0;
			}
			for (unsigned int i=isXmin ; i<nl_d ; i++) {
				(*BtAM)(i,j)= -(*BtAM)(i,j+1);
			}
			for (unsigned int i=isXmin ; i<nl_p ; i++) {
				(*BlAM)(i,j)= (*BlAM)(i,j+1);
				//(*BlAM)(i,0)+= -(*BlAM)(i,1)+(*BlAM_old)(i,1)-4*dt_ov_dr*(*EtAM)(i,1);
			}
		}

		else if (imode==1){
			//MF
			for (unsigned int i=isXmin ; i<nl_p  ; i++) {
				(*BlAM)(i,j)= -(*BlAM)(i,j+1);
                //if (std::abs((*BlAM)(i,j))>1.){
                //MESSAGE("BlAMA");                
                //MESSAGE(i);
                //MESSAGE(j);    
                //MESSAGE((*BlAM)(i,j));
                //}
			}

			for (unsigned int i=isXmin ; i<nl_d-1 ; i++) {
				(*BrAM)(i,j)+=  Icpx*dt_ov_dr*(*ElAM)(i,j+1)
				+			dt_ov_dl*((*EtAM)(i,j)-(*EtAM)(i-1,j));
                //if (std::abs((*BrAM)(i,j))>1.){
                //MESSAGE("BrAMA");                
                //MESSAGE(i);
                //MESSAGE(j);    
                //MESSAGE((*BrAM)(i,j));
                //}
			}
			for (unsigned int i=isXmin ; i<nl_d ; i++) {
				//(*BtAM)(i,0)+= -dt_ov_dl*((*ErAM)(i+1,0)-(*ErAM)(i,0)+(*ErAM)(i+1,1)-(*ErAM)(i,1))
				//+				2*dt_ov_dr*(*ElAM)(i+1,1) - (*BtAM_old)(i,1)+ (*BtAM)(i,1);
				(*BtAM)(i,j)= -2.*Icpx*(*BrAM)(i,j)-(*BtAM)(i,j+1);
                //if (std::abs((*BtAM)(i,j))>1.){
                //MESSAGE("BtAMA");                
                //MESSAGE(i);
                //MESSAGE(j);    
                //MESSAGE((*BtAM)(i,j));
                //}
			}	

		}
		else {
			for (unsigned int  i=isXmin ; i<nl_d; i++) {
				(*BlAM)(i,j)= -(*BlAM)(i,j+1);
			}
			for (unsigned int i=isXmin ; i<nl_p; i++) {
				(*BrAM)(i,j)= 0;
			}
			for (unsigned int  i=isXmin ; i<nl_d ; i++) {
				(*BtAM)(i,j)= - (*BtAM)(i,j+1);
			}	

		}
	}    
	}
}

