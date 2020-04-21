#include "EnvelopeBCAM_Axis.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "Field2D.h"
#include "cField2D.h"
#include "Tools.h"
#include "LaserEnvelope.h"
#include "ElectroMagn.h"

using namespace std;


EnvelopeBCAM_Axis::EnvelopeBCAM_Axis( Params &params, Patch *patch, unsigned int _min_max )
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
void EnvelopeBCAM_Axis::apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    // Static cast of the field
    cField2D *A2Dcyl        = static_cast<cField2D *>( envelope->A_ );  // the envelope at timestep n
    
    Field2D  *Phi2Dcyl      = static_cast<Field2D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n
    
    //Field2D *Env_Aabs2Dcyl  = static_cast<Field2D *>( EMfields->Env_A_abs_ ); // absolute value of the envelope |A|

    Field2D *Env_Eabs2Dcyl  = static_cast<Field2D *>( EMfields->Env_E_abs_ ); // absolute value of the envelope of the transverse electric field of the  laser |E|
 
    Field2D *Env_Exabs2Dcyl  = static_cast<Field2D *>( EMfields->Env_Ex_abs_ ); // absolute value of the envelope of the longitudinal electric field of the  laser |Ex|
    
    double ellipticity_factor = envelope->ellipticity_factor;
// APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
  

    if( min_max == 2 && patch->isYmin() ) { // j_p = 2 corresponds to r=0
    
        // zero radial derivative on axis 
        //unsigned int j=2;
        for( unsigned int i=0; i<nx_p; i++ ) {
           ( *A2Dcyl )         ( i, 1 ) = (*A2Dcyl)(i,3);
           ( *Phi2Dcyl )       ( i, 1 ) = ellipticity_factor*std::abs((*A2Dcyl)  (i,3)) * std::abs((*A2Dcyl)  (i,3)) * 0.5;
           // ( *Env_Aabs2Dcyl )  ( i, 1 ) = (*Env_Aabs2Dcyl)(i,3);
           // ( *Env_Eabs2Dcyl )  ( i, 1 ) = (*Env_Eabs2Dcyl)(i,3);
           ( *A2Dcyl )         ( i, 0 ) =  (*A2Dcyl)(i,4);
           ( *Phi2Dcyl )       ( i, 0 ) = ellipticity_factor*std::abs((*A2Dcyl)  (i,4)) * std::abs((*A2Dcyl)  (i,4)) * 0.5;
           // ( *Env_Aabs2Dcyl )  ( i, 0 ) = (*Env_Aabs2Dcyl)(i,4);
           // ( *Env_Eabs2Dcyl )  ( i, 0 ) = (*Env_Eabs2Dcyl)(i,4);

            
        }//i
        
    } 
    
    
    
    
    
    
}

