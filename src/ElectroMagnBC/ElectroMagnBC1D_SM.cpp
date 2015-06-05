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

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply_xmin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
    if ( smpi->isWestern() ) {
        
        Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
        
        // Laser temporal profile
        double byL=0.0, bzL=0.0;
        
        for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
            
            if (laser_[ilaser]->laser_struct.boxSide == "west") {
                // omega*time
                double wt = laser_[ilaser]->omega0_ * time_dual
                *         ( 1.0 + laser_[ilaser]->tchirp_*laser_[ilaser]->omega0_*time_dual );
                // Incident field (left boundary)
                byL += laser_[ilaser]->a0_delta_y_ * sin(wt) * laser_[ilaser]->time_profile(time_dual);
                bzL += laser_[ilaser]->a0_delta_z_ * cos(wt) * laser_[ilaser]->time_profile(time_dual);
            } else if (laser_[ilaser]->laser_struct.boxSide == "east") {
                // nothing done (treating only west bnd here)
            } else {
                ERROR("Angle not allowed for 1D/2D laser pulse " << ilaser);
            }
            
        }//ilaser
        
        // Apply Silver-Mueller EM boundary condition at x=xmin
        (*By1D)(0) =  Alpha_SM*(*Ez1D)(0) + Beta_SM*(*By1D)(1) + Gamma_SM*byL;
        (*Bz1D)(0) = -Alpha_SM*(*Ey1D)(0) + Beta_SM*(*Bz1D)(1) + Gamma_SM*bzL;
        
    }//if Western

}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply_xmax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
    
    if ( smpi->isEastern() ) {
        Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
    
        double byR=0.0, bzR=0.0;
    
        for (unsigned int ilaser=0; ilaser< laser_.size(); ilaser++) {
        
            if (laser_[ilaser]->laser_struct.boxSide == "west") {
                // nothing done (treating only east bnd here)
            } else if (laser_[ilaser]->laser_struct.boxSide == "east") {
                // omega*time
                double wt = laser_[ilaser]->omega0_ * time_dual
                * ( 1.0 + laser_[ilaser]->tchirp_*laser_[ilaser]->omega0_*time_dual );
                // Incident field (right boundary)
                byR += laser_[ilaser]->a0_delta_y_ * sin(wt) * laser_[ilaser]->time_profile(time_dual);
                bzR += laser_[ilaser]->a0_delta_z_ * cos(wt) * laser_[ilaser]->time_profile(time_dual);
            } else {
                ERROR("Angle not allowed for 1D/2D laser pulse " << ilaser);
            }
        
        }//ilaser
    
        // Silver-Mueller boundary conditions (right)
        (*By1D)(nx_d-1) = -Alpha_SM*(*Ez1D)(nx_p-1) + Beta_SM*(*By1D)(nx_d-2) + Gamma_SM*byR;
        (*Bz1D)(nx_d-1) =  Alpha_SM*(*Ey1D)(nx_p-1) + Beta_SM*(*Bz1D)(nx_d-2) + Gamma_SM*bzR;
        
    }//if Eastern
    
}

void ElectroMagnBC1D_SM::apply_ymin(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
}
void ElectroMagnBC1D_SM::apply_ymax(ElectroMagn* EMfields, double time_dual, SmileiMPI* smpi)
{
}
