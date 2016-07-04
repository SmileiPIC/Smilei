#include "ElectroMagnBC1D_SM.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Tools.h"
#include "Laser.h"

using namespace std;

ElectroMagnBC1D_SM::ElectroMagnBC1D_SM( Params &params, Patch* patch )
  : ElectroMagnBC( params, patch )
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
    
    
    Bz_xvalmin = 0.;
    Bz_xvalmax = 0.;
    By_xvalmin = 0.;
    By_xvalmax = 0.;
}

ElectroMagnBC1D_SM::~ElectroMagnBC1D_SM()
{
}

void ElectroMagnBC1D_SM::save_fields_BC1D(Field* my_field) {
    Field1D* field1D=static_cast<Field1D*>(my_field);
    // Bx^(p) is not saved as it is defined on the primal grid and thus can be computed
    // we save only the field By and Bz that are computed on the dual grid
    
    if (field1D->name=="By"){
        By_xvalmin=(*field1D)(0);
        By_xvalmax=(*field1D)(field1D->dims()[0]-1);
    }
    
    if (field1D->name=="Bz"){
        Bz_xvalmin=(*field1D)(0);
        Bz_xvalmax=(*field1D)(field1D->dims()[0]-1);
    }

}



// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    if ( patch->isWestern() ) {
        
        //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
        
        // Lasers
        double byL=0.0, bzL=0.0;
        vector<double> pos(1);
        pos[0] = 0.;
        for (int ilaser=0; ilaser<vecLaser.size(); ilaser++) {
            byL += vecLaser[ilaser]->getAmplitude0(pos, time_dual, 0);
            bzL += vecLaser[ilaser]->getAmplitude1(pos, time_dual, 0);
        }
        
        // Apply Silver-Mueller EM boundary condition at x=xmin
        (*By1D)(0) =  Alpha_SM*(*Ez1D)(0) + Beta_SM*((*By1D)(1)-By_xvalmin) + Gamma_SM*byL+By_xvalmin;
        (*Bz1D)(0) = -Alpha_SM*(*Ey1D)(0) + Beta_SM*((*Bz1D)(1)-Bz_xvalmin) + Gamma_SM*bzL+Bz_xvalmin;
        
    }//if Western

}
// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_SM::apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    
    if ( patch->isEastern() ) {
        //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
    
        // Lasers
        double byR=0.0, bzR=0.0;
        vector<double> pos(1);
        pos[0] = 0.;
        for (unsigned int ilaser=0; ilaser< vecLaser.size(); ilaser++) {
            byR += vecLaser[ilaser]->getAmplitude0(pos, time_dual, 0);
            bzR += vecLaser[ilaser]->getAmplitude1(pos, time_dual, 0);
        }
    
        // Silver-Mueller boundary conditions (right)
        (*By1D)(nx_d-1) = -Alpha_SM*(*Ez1D)(nx_p-1)+ Beta_SM*((*By1D)(nx_d-2)-By_xvalmax) + Gamma_SM*byR+By_xvalmax;
        (*Bz1D)(nx_d-1) =  Alpha_SM*(*Ey1D)(nx_p-1)+ Beta_SM*((*Bz1D)(nx_d-2)-Bz_xvalmax) + Gamma_SM*bzR+Bz_xvalmax;
    }//if Eastern
    
}

void ElectroMagnBC1D_SM::apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}
void ElectroMagnBC1D_SM::apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}

void ElectroMagnBC1D_SM::apply_zmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}
void ElectroMagnBC1D_SM::apply_zmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}
