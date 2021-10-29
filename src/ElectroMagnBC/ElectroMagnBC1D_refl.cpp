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
ElectroMagnBC1D_refl::ElectroMagnBC1D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : ElectroMagnBC1D( params, patch, i_boundary )
{

    // oversize
    oversize_ = params.oversize[0];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Reflective Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void ElectroMagnBC1D_refl::apply( ElectroMagn *EMfields, double time_dual, Patch *patch )
{
    if( i_boundary_ == 0 ) {
        if( patch->isXmin() ) {
        
            // Application over the full-ghost cell
            //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
            //Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
            //Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
            //Field1D* Bx1D   = static_cast<Field1D*>(EMfields->Bx_);
            Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
            Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
            
            // force constant magnetic fields in the ghost cells
            for( unsigned int i=oversize_; i>0; i-- ) {
                //(*Bx1D)(i-1) = (*Bx1D)(i);
                ( *By1D )( i-1 ) = ( *By1D )( i );
                ( *Bz1D )( i-1 ) = ( *Bz1D )( i );
            }
            
            //        // force 0 electric fields in the ghost cells
            //        for (unsigned int i=0; i<oversize_+1; i++) {
            //            (*Ex1D)(i) = 0.0;
            //        }
            //        for (unsigned int i=0; i<oversize_; i++) {
            //            (*Ey1D)(i) = 0.0;
            //            (*Ez1D)(i) = 0.0;
            //        }
            
            
            /* DEFINITION BY NICO
            
             // The other fields are already defined everywhere, so no need for boundary conditions.
             Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
             Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
            
             // normal derivative of tangential B is zero.
             // By and Bz just outside equal By and Bz just inside.
             (*By1D)(0) = (*By1D)(1);
             (*Bz1D)(0) = (*Bz1D)(1);
             */
            
        }//if Xmin
    } else {
        if( patch->isXmax() ) {
        
            // application of Bcs over the full ghost cells
            //Field1D* Ex1D   = static_cast<Field1D*>(EMfields->Ex_);
            //Field1D* Ey1D   = static_cast<Field1D*>(EMfields->Ey_);
            //Field1D* Ez1D   = static_cast<Field1D*>(EMfields->Ez_);
            //Field1D* Bx1D   = static_cast<Field1D*>(EMfields->Bx_);
            Field1D *By1D   = static_cast<Field1D *>( EMfields->By_ );
            Field1D *Bz1D   = static_cast<Field1D *>( EMfields->Bz_ );
            
            // force constant magnetic fields in the ghost cells
            //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++)
            //            (*Bx1D)(i) = (*Bx1D)(i-1);
            for( unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++ ) {
                ( *By1D )( i ) = ( *By1D )( i-1 );
                ( *Bz1D )( i ) = ( *Bz1D )( i-1 );
            }
            
            //        // force 0 electric fields in the ghost cells
            //        for (unsigned int i=n_d[0]-oversize_; i<n_d[0]; i++)
            //            (*Ex1D)(i) = 0.0;
            //        for (unsigned int i=n_p[0]-oversize_; i<n_p[0]; i++) {
            //            (*Ey1D)(i) = 0.0;
            //            (*Ez1D)(i) = 0.0;
            //        }
            
            /* DEFINITION BY NICO
             Field1D* By1D   = static_cast<Field1D*>(EMfields->By_);
             Field1D* Bz1D   = static_cast<Field1D*>(EMfields->Bz_);
            
             // normal derivative of tangential B is zero.
             // By and Bz just outside equal By and Bz just inside.
             (*By1D)(n_d[0]-1) = (*By1D)(n_d[0]-2);
             (*Bz1D)(n_d[0]-1) = (*Bz1D)(n_d[0]-2);
             */
            
        }//if Xmax
        
    }
}

