#include <cstdlib>
#include <iostream>
#include <string>

#include "EnvelopeBC.h"

#include "Params.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;

// Constructor for EnvelopeBC
EnvelopeBC::EnvelopeBC( Params &params, Patch *, unsigned int i_boundary ) :
    i_boundary_( i_boundary )
{

    // time step
    dt = params.timestep;

    std::vector<unsigned int> n_space(params.patch_size_);
    if (params.multiple_decomposition)
        n_space = params.region_size_;

    std::vector<unsigned int> oversize(params.oversize);
    if (params.multiple_decomposition)
        oversize = params.region_oversize;

    // number of nodes of the primal and dual grid in the x-direction
    nx_p = n_space[0]+1+2*oversize[0];
    nl_p = n_space[0]+1+2*oversize[0];
    // number of nodes of the primal and dual grid in the y-direction
    if (params.nDim_field>1){
        ny_p = n_space[1]+1+2*oversize[1];
        nr_p = n_space[1]+1+2*oversize[1];
    }
    // number of nodes of the primal and dual grid in the z-direction
    if (params.nDim_field>2){
        nz_p = n_space[2]+1+2*oversize[2];
    }

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the x-direction)
    dx       = params.cell_length[0];
    dl       = params.cell_length[0];
    dt_ov_dx = dt/dx;
    dx_ov_dt = 1.0/dt_ov_dx;

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the y-direction)
    if (params.nDim_field>1){
        dy       = params.cell_length[1];
        dr       = params.cell_length[1];
        dt_ov_dy = dt/dy;
        dy_ov_dt = 1.0/dt_ov_dy;
    }

    // spatial-step and ratios time-step by spatial-step & spatial-step by time-step (in the z-direction)
    if (params.nDim_field>2){
        dz       = params.cell_length[2];
        dt_ov_dz = dt/dz;
        dz_ov_dt = 1.0/dt_ov_dz;
    }

    pml_solver_envelope_ = NULL;
}

// Destructor for EnvelopeBC
EnvelopeBC::~EnvelopeBC()
{
}
