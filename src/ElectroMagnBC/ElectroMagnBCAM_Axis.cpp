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

ElectroMagnBCAM_Axis::ElectroMagnBCAM_Axis( Params &params, Patch *patch, unsigned int _min_max )
    : ElectroMagnBCAM( params, patch, _min_max )
{
    // conversion factor from degree to radian
    conv_deg2rad = M_PI/180.0;
    
    //Number of modes
    Nmode= params.nmodes;
}


void ElectroMagnBCAM_Axis::save_fields( Field *my_field, Patch *patch )
{
//    ERROR( "Impossible" );
}

void ElectroMagnBCAM_Axis::disableExternalFields()
{
    ERROR( "Impossible" );
}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBCAM_Axis::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    return; // For the moment, boundary conditions on axis are handled directly in the solvers.
}

