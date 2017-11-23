#include "ElectroMagnBCRZ_Axis.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"
#include "Tools.h"

using namespace std;

ElectroMagnBCRZ_Axis::ElectroMagnBCRZ_Axis( Params &params, Patch* patch, unsigned int _min_max )
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
    
    if (min_max == 2 && patch->isYmin() ) {
    }
    
    #ifdef _TODO_RZ
    #endif
}


void ElectroMagnBCRZ_Axis::save_fields(Field* my_field, Patch* patch)
{
    ERROR( "Impossible" );
}

void ElectroMagnBCRZ_Axis::disableExternalFields()
{
    ERROR( "Impossible" );
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCRZ_Axis::apply(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    // Loop on imode 
    int imode = 0;

    // Static cast of the fields
    cField2D* El2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->El_[imode];
    cField2D* Er2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->Er_[imode];
    cField2D* Et2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->Et_[imode];
    cField2D* Bl2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bl_[imode];
    cField2D* Br2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->Br_[imode];
    cField2D* Bt2D = (static_cast<ElectroMagn3DRZ*>(EMfields))->Bt_[imode];
 
    if (min_max == 2 && patch->isYmin() ) {
    }
    
}
