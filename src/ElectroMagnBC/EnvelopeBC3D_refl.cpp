#include "EnvelopeBC3D_refl.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Params.h"
#include "Patch.h"
#include "Field3D.h"
#include "cField3D.h"
#include "Tools.h"
#include "LaserEnvelope.h"
#include "ElectroMagn.h"

using namespace std;


EnvelopeBC3D_refl::EnvelopeBC3D_refl( Params &params, Patch *patch, unsigned int i_boundary )
    : EnvelopeBC( params, patch, i_boundary )
{
    // oversize
    oversize_ = params.oversize[0];
    
    // number of nodes of the primal and dual grid in the x-direction
    nx_p = params.n_space[0]+1+2*params.oversize[0];
    nx_d = params.n_space[0]+2+2*params.oversize[0];
    
    // number of nodes of the primal and dual grid in the y-direction
    ny_p = params.n_space[1]+1+2*params.oversize[1];
    ny_d = params.n_space[1]+2+2*params.oversize[1];
    
    // number of nodes of the primal and dual grid in the z-direction
    nz_p = params.n_space[2]+1+2*params.oversize[2];
    nz_d = params.n_space[2]+2+2*params.oversize[2];
    
}

// ---------------------------------------------------------------------------------------------------------------------
// Apply Boundary Conditions
// ---------------------------------------------------------------------------------------------------------------------
void EnvelopeBC3D_refl::apply( LaserEnvelope *envelope, ElectroMagn *EMfields, double time_dual, Patch *patch )
{

    // Static cast of the fields immediately after the envelope equation sover
    cField3D *A3D        = static_cast<cField3D *>( envelope->A_ );  // the envelope at timestep n+1
    Field3D  *Phi3D      = static_cast<Field3D *>( envelope->Phi_ ); // the ponderomotive potential Phi=|A|^2/2 at timestep n+1
    
    // APPLICATION OF BCs OVER THE FULL GHOST CELL REGION
    
    if( i_boundary_ == 0 && patch->isXmin() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=oversize_; i>0; i-- ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    ( *A3D )( i-1, j, k ) = 0. ; // (*A3D)(i,j,k);
                    ( *Phi3D )( i-1, j, k ) = 0. ; // std::abs((*A3D)  (i,j,k)) * std::abs((*A3D)  (i,j,k)) * 0.5;
                }//k
            }//j
        }//i
        
    } else if( i_boundary_ == 1 && patch->isXmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=nx_p-oversize_; i<nx_p; i++ ) {
            for( unsigned int j=0 ; j<ny_p ; j++ ) {
                for( unsigned int k=0 ; k<nz_p ; k++ ) {
                    ( *A3D )( i, j, k ) = 0. ; //(*A3D)(i-1,j,k);
                    ( *Phi3D )( i, j, k ) = 0. ; //std::abs((*A3D)  (i-1,j,k)) * std::abs((*A3D)  (i-1,j,k)) * 0.5;
                }//k
            }//j
        }//i
        
    } else if( i_boundary_ == 2 && patch->isYmin() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=0; i<nx_p; i++ ) {
            for( unsigned int j=oversize_ ; j>0 ; j-- ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *A3D )( i, j-1, k ) = 0. ; // (*A3D)(i,j,k);
                    ( *Phi3D )( i, j-1, k ) = 0. ; // std::abs((*A3D)  (i,j,k)) * std::abs((*A3D)  (i,j,k)) * 0.5;
                }//k
            }//j
        }//i
        
    } else if( i_boundary_ == 3 && patch->isYmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=0; i<nx_p; i++ ) {
            for( unsigned int j=ny_p-oversize_; j<ny_p ; j++ ) {
                for( unsigned int k=0; k<nz_p; k++ ) {
                    ( *A3D )( i, j, k ) = 0. ; // (*A3D)(i,j-1,k);
                    ( *Phi3D )( i, j, k ) = 0. ; // std::abs((*A3D)  (i,j-1,k)) * std::abs((*A3D)  (i,j-1,k)) * 0.5;
                }//k
            }//j
        }//i
        
    } else if( i_boundary_ == 4 && patch->isZmin() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=0; i<nx_p; i++ ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=oversize_; k>0; k-- ) {
                    ( *A3D )( i, j, k-1 ) = 0. ; // (*A3D)(i,j,k);
                    ( *Phi3D )( i, j, k-1 ) = 0. ; // std::abs((*A3D)  (i,j,k)) * std::abs((*A3D)  (i,j,k)) * 0.5;
                }//k
            }//j
        }//i
        
    } else if( i_boundary_ == 5 && patch->isZmax() ) {
    
        // FORCE CONSTANT ENVELOPE FIELD ON BORDER
        
        for( unsigned int i=0; i<nx_p; i++ ) {
            for( unsigned int j=0; j<ny_p; j++ ) {
                for( unsigned int k=nz_p-oversize_; k<nz_p; k++ ) {
                    ( *A3D )( i, j, k ) = 0. ; // (*A3D)(i,j,k-1);
                    ( *Phi3D )( i, j, k ) = 0. ; // std::abs((*A3D)  (i,j,k-1)) * std::abs((*A3D)  (i,j,k-1)) * 0.5;
                }//z
            }//j
        }//i
        
    }
    
    
    
    
}

