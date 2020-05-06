
#include "ElectroMagnBCAM.h"

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "cField2D.h"


ElectroMagnBCAM::ElectroMagnBCAM( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBC( params, patch, _min_max )
{
    std::vector<unsigned int> n_space(params.n_space);
    std::vector<unsigned int> oversize(params.oversize);
    if (params.uncoupled_grids) {
        n_space = params.n_space_region;
        oversize = params.region_oversize;
    }

    // number of nodes of the primal and dual grid in the l-direction
    nl_p = n_space[0]+1+2*oversize[0];
    nl_d = nl_p+1-params.is_pxr;
    // number of nodes of the primal and dual grid in the r-direction
    nr_p = n_space[1]+1+2*oversize[1];
    nr_d = nr_p+1-params.is_pxr;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the l-direction)
    dl       = params.cell_length[0];
    dt_ov_dl = dt/dl;
    dl_ov_dt = 1.0/dt_ov_dl;
    
    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    dr       = params.cell_length[1];
    dt_ov_dr = dt/dr;
    dr_ov_dt = 1.0/dt_ov_dr;
    
}

ElectroMagnBCAM::~ElectroMagnBCAM()
{
}


void ElectroMagnBCAM::applyBConEdges( ElectroMagn *EMfields, Patch *patch )
{
}
