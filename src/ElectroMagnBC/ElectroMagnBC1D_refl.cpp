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
        
        Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* Bx1D   = static_cast<Field1D*>(EMfields->Bx_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
        
        // force constant magnetic fields in the ghost cells
        for (unsigned int i=oversize_; i>0; i-- ) {
            (*Bx1D)(i-1) = (*Bx1D)(i);
            (*By1D)(i-1) = (*By1D)(i);
            (*Bz1D)(i-1) = (*Bz1D)(i);
        }
        
        // force 0 electric fields in the ghost cells
        for (unsigned int i=0; i<oversize_+1; i++) {
            (*Ex1D)(i) = 0.0;
        }
        for (unsigned int i=0; i<oversize_; i++) {
            (*Ey1D)(i) = 0.0;
            (*Ez1D)(i) = 0.0;
        }
        
    }//if Western

}


// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply_xmax(ElectroMagn* EMfields, double time_dual, Patch* patch)
{
    
    if ( patch->isEastern() ) {
        
        Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
        Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
        Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
        Field1D* Bx1D   = static_cast<Field1D*>(EMfields->Bx_);
        Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
        Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
        
        // force constant magnetic fields in the ghost cells
        for (unsigned int i=nx_p-oversize_; i<nx_p; i++)
            (*Bx1D)(i) = (*Bx1D)(i-1);
        for (unsigned int i=nx_d-oversize_; i<nx_d; i++) {
            (*By1D)(i) = (*By1D)(i-1);
            (*Bz1D)(i) = (*Bz1D)(i-1);
        }
        
        // force 0 electric fields in the ghost cells
        for (unsigned int i=nx_d-oversize_; i<nx_d; i++)
            (*Ex1D)(i) = 0.0;
        for (unsigned int i=nx_p-oversize_; i<nx_p; i++) {
            (*Ey1D)(i) = 0.0;
            (*Ez1D)(i) = 0.0;
        }
        
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
