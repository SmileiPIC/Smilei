#include "EnvelopeBCAM_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "Field2D.h"
#include "cField2D.h"
#include "Tools.h"
#include "LaserEnvelope.h"

using namespace std;


EnvelopeBCAM_refl::EnvelopeBCAM_refl( Params &params, Patch *patch, unsigned int _min_max )
    : EnvelopeBC( params, patch, _min_max )
{
    // oversize
    oversize_ = params.oversize[0];
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = params.n_space[0]+2+2*params.oversize[0];
    
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]+1+2*params.oversize[1];
    ny_d = params.n_space[1]+2+2*params.oversize[1];
    
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void EnvelopeBCAM_refl::apply( LaserEnvelope *envelope, double time_dual, Patch *patch )
{

    // Static cast of the field
    cField2D *A2Dcyl        = static_cast<cField2D *>( envelope->A_ );  // the envelope at timestep n
    
    Field2D  *Phi2Dcyl      = static_cast<Field2D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n
    
    
    // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
    
    if( min_max == 0 && patch->isXmin() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=oversize_; i>0; i-- ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                ( *A2Dcyl )( i-1, j ) = 0. ; // (*A2D)(i,j,k);
                ( *Phi2Dcyl )( i-1, j ) = 0. ; // std::abs((*A2D)  (i,j)) * std::abs((*A2D)  (i,j)) * 0.5;
            }//j
        }//i
        
    } else if( min_max == 1 && patch->isXmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=nx_p-oversize_; i<nx_p; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                ( *A2Dcyl )( i, j ) = 0. ; //(*A2D)(i-1,j);
                ( *Phi2Dcyl )( i, j ) = 0. ; //std::abs((*A2D)  (i-1,j)) * std::abs((*A2D)  (i-1,j)) * 0.5;
            }//j
        }//i
        
    } else if( min_max == 3 && patch->isYmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=0; i<nx_p; i++ ) {
            for( unsigned int j=ny_p-oversize_; j<ny_p ; j++ ) {
                ( *A2Dcyl )( i, j ) = 0. ; // (*A2D)(i,j-1);
                ( *Phi2Dcyl )( i, j ) = 0. ; // std::abs((*A2D)  (i,j-1)) * std::abs((*A2D)  (i,j-1)) * 0.5;
            }//j
        }//i
        
    }
    
    
    
    
    
    
}

