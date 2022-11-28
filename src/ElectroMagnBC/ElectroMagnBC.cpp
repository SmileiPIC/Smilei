
#include <cstdlib>
#include <iostream>
#include <string>

#include "ElectroMagnBC.h"

#include "Params.h"
#include "Laser.h"
#include "Tools.h"
#include "Patch.h"

using namespace std;

// Constructor for ElectromagnBC
ElectroMagnBC::ElectroMagnBC( Params &params, Patch *patch, unsigned int i_boundary ) :
    i_boundary_( i_boundary )
{
    vecLaser.resize( 0 );
    
    // time step
    dt = params.timestep;
    
    n_p.resize( params.nDim_field );
    n_d.resize( params.nDim_field );
    d.resize( params.nDim_field );
    dt_ov_d.resize( params.nDim_field );
    for( unsigned int i=0; i<params.nDim_field; i++ ) {
        // number of nodes of the primal and dual grid
        n_p[i] = patch->size_[i] + 1 + 2*patch->oversize[i];
        n_d[i] = n_p[i] + 1 - params.is_pxr;
        
        // spatial-step and ratios time-step by spatial-step
        d[i] = params.cell_length[i];
        dt_ov_d[i] = dt / d[i];
    }

    pml_solver_ = NULL;
}

// Destructor for ElectromagnBC
ElectroMagnBC::~ElectroMagnBC()
{
    for( unsigned int i=0; i< vecLaser.size(); i++ ) {
        delete vecLaser[i];
    }
    vecLaser.clear();
}


// Disable all lasers when using moving window
void ElectroMagnBC::laserDisabled()
{
    //for (unsigned int i=0; i< vecLaser.size(); i++) {
    //    vecLaser[i]->disable();
    //}
    for( unsigned int i=0; i< vecLaser.size(); i++ ) {
        delete vecLaser[i];
    }
    vecLaser.resize( 0 );
}

