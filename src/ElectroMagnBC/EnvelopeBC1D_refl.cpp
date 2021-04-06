#include "EnvelopeBC1D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "Field1D.h"
#include "cField1D.h"
#include "Tools.h"
#include "LaserEnvelope.h"
#include "ElectroMagn.h"

using namespace std;


EnvelopeBC1D_refl::EnvelopeBC1D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : EnvelopeBC( params, patch, i_boundary )
{
    // oversize
    oversize_ = params.oversize[0];
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = params.n_space[0]+2+2*params.oversize[0];
    
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void EnvelopeBC1D_refl::apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    // Static cast of the field
    cField1D *A1D        = static_cast<cField1D *>( envelope->A_ );  // the envelope at timestep n
    Field1D  *Phi1D      = static_cast<Field1D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n
    
    // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
    
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=oversize_; i>0; i-- ) {
            ( *A1D )( i-1 ) = 0. ;  // (*A1D)(i,j,k);
            ( *Phi1D )( i-1 ) = 0. ; // std::abs((*A1D)  (i)) * std::abs((*A1D)  (i)) * 0.5;
        }//i
        
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=nx_p-oversize_; i<nx_p; i++ ) {
            ( *A1D )( i ) = 0. ;  //(*A1D)(i-1,j);
            ( *Phi1D )( i ) = 0. ; //std::abs((*A1D)  (i-1)) * std::abs((*A1D)  (i-1)) * 0.5;
        }//i
        
    }
    
    
}

