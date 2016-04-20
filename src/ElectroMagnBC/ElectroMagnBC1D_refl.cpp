#include "ElectroMagnBC1D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "ElectroMagn.h"
#include "Field1D.h"
#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Reflective boundary conditions
// The ghost cells behave as a Perfect Condutor:
//      - electric fields are put to zero in the ghost cells
//      - magnetic fields are constant in the ghost cells
// ---------------------------------------------------------------------------------------------------------------------
ElectroMagnBC1D_refl::ElectroMagnBC1D_refl( Params &params, Patch* patch )
  : ElectroMagnBC( params, patch )
{

    // oversize
    oversize_ = params.oversize[0];

    // number of nodes of the primal-grid
    nx_p = params.n_space[0]+1 + 2*params.oversize[0];

    // number of nodes of the dual-grid
    nx_d = params.n_space[0]+2 + 2*params.oversize[0];

}

ElectroMagnBC1D_refl::~ElectroMagnBC1D_refl()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply_xmin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{

    if ( patch->isWestern() ) {
        // The other fields are already defined everywhere, so no need for boundary conditions.
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);

        // normal derivative of tangential B is zero.
        // By and Bz just outside equal By and Bz just inside.
        (*By1D)(0) = (*By1D)(1);
        (*Bz1D)(0) = (*Bz1D)(1);
    }//if Western

}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{

    if ( patch->isEastern() ) {
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);

        // normal derivative of tangential B is zero.
        // By and Bz just outside equal By and Bz just inside.
        (*By1D)(nx_d-1) = (*By1D)(nx_d-2);
        (*Bz1D)(nx_d-1) = (*Bz1D)(nx_d-2);

    }//if Eastern

}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply_ymin(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}
void ElectroMagnBC1D_refl::apply_ymax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
}
