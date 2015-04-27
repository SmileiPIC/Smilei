#include "ElectroMagnBC1D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Tools.h"

using namespace std;

ElectroMagnBC1D_SM::ElectroMagnBC1D_SM( PicParams &params, LaserParams &laser_params )
    : ElectroMagnBC( params, laser_params )
{
    // number of nodes of the primal-grid
    nx_p = params.n_space[0]+1 + 2*params.oversize[0];
    // number of nodes of the dual-grid
    nx_d = params.n_space[0]+2 + 2*params.oversize[0];

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step
    dx       = params.cell_length[0];
    dt_ov_dx = params.timestep/params.cell_length[0];
    dx_ov_dt = 1.0/dt_ov_dx;

    // Parameters for the Silver-Mueller boundary conditions
    Alpha_SM = 2./(1.+dt_ov_dx);
    Beta_SM  = (dt_ov_dx-1.)/(1.+dt_ov_dx);
    Gamma_SM = 4./(1.+dt_ov_dx);

}

ElectroMagnBC1D_SM::~ElectroMagnBC1D_SM()
{
}

void ElectroMagnBC1D_SM::save_fields_value_for_BC(Field* my_field) {
    Field1D* field1D=static_cast<Field1D*>(my_field);
    xvalmin=(*field1D)(0);
    xvalmax=(*field1D)(field1D->dims()[0]-1);
}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
//    Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
    Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
    Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
    Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
    Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);

    // --------------------------------------------------
    // Laser temporal profile
    // --------------------------------------------------
    double byL=0.0, bzL=0.0, byR=0.0, bzR=0.0;

    for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {

        if (laser_[ilaser]->laser_struct.boxSide == "west") {
            // Incident field (left boundary)
            byL += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
            bzL += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
        } else if (laser_[ilaser]->laser_struct.boxSide == "east") {
            // Incident field (right boundary)
            byR += laser_[ilaser]->a0_delta_y_ * sin(time_dual) * laser_[ilaser]->time_profile(time_dual);
            bzR += laser_[ilaser]->a0_delta_z_ * cos(time_dual) * laser_[ilaser]->time_profile(time_dual);
        } else {
            ERROR("Angle not allowed for 1D/2D laser pulse " << ilaser);
        }

    }//ilaser

    // ----------------------------
    // Apply EM boundary conditions
    // In case of External field, we apply the boundary conditions on the field created by the plasma -> Bz-Bzext
    // where Bzext is the external field defined as initial conditon 
    // ----------------------------
    if ( smpi->isWestern() ) {
        // Silver-Mueller boundary conditions (left)
        (*By1D)(0) =  Alpha_SM*((*Ez1D)(0)-Ez1D->ExtFieldAt_xmin) + Beta_SM*((*By1D)(1)-By1D->ExtFieldAt_xmin) + Gamma_SM*byL;
        (*Bz1D)(0) = -Alpha_SM*((*Ey1D)(0)-Ey1D->ExtFieldAt_xmin) + Beta_SM*((*Bz1D)(1)-Bz1D->ExtFieldAt_xmin) + Gamma_SM*bzL;
    }//if Western
    if ( smpi->isEastern() ) {
        // Silver-Mueller boundary conditions (right)
        (*By1D)(nx_d-1) = -Alpha_SM*((*Ez1D)(nx_p-1)-Ez1D->ExtFieldAt_xmax) + Beta_SM*((*By1D)(nx_d-2)-By1D->ExtFieldAt_xmax) + Gamma_SM*byR;
        (*Bz1D)(nx_d-1) =  Alpha_SM*((*Ey1D)(nx_p-1)-Ey1D->ExtFieldAt_xmax) + Beta_SM*((*Bz1D)(nx_d-2)-Bz1D->ExtFieldAt_xmax) + Gamma_SM*bzR;
    }//if Eastern

}
