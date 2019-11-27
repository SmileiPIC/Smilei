#include "Pusher.h"
#include "Params.h"
#include "Species.h"

Pusher::Pusher( Params &params, Species *species ) :
    min_loc_vec( species->min_loc_vec ),
    vecto( params.vectorization_mode=="on" || params.vectorization_mode=="adaptive_mixed_sort" || params.vectorization_mode=="adaptive" || params.cell_sorting )
{
    for( unsigned int ipos=0; ipos < params.nDim_particle ; ipos++ ) {
        dx_inv_[ipos] = species->dx_inv_[ipos];
    }
    
    nspace[0] = 0;
    nspace[1] = params.n_space[1]+1;
    nspace[2] = params.n_space[2]+1;
    
    mass_          = species->mass_;
    if( mass_ > 0. ) {
        one_over_mass_ = 1.0/mass_;
    } else {
        one_over_mass_ = 0.;
    }
    dt             = params.timestep;
    dts2           = params.timestep/2.;
    dts4           = params.timestep/4.;
    
    nDim_          = params.nDim_particle;
    
}

Pusher::~Pusher()
{
}
