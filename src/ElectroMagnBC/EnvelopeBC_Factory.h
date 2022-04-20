#ifndef ENVELOPEBC_FACTORY_H
#define ENVELOPEBC_FACTORY_H

#include "EnvelopeBC.h"
#include "EnvelopeBCAM_Axis.h"
#include "EnvelopeBCAM_refl.h"
#include "EnvelopeBCAM_PML.h"
#include "EnvelopeBC3D_refl.h"
#include "EnvelopeBC3D_PML.h"
#include "EnvelopeBC2D_refl.h"
#include "EnvelopeBC2D_PML.h"
#include "EnvelopeBC1D_refl.h"

#include "Params.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Constructor for the Envelope Boundary conditions factory
//  --------------------------------------------------------------------------------------------------------------------
class EnvelopeBC_Factory
{

public:

    static std::vector<EnvelopeBC *> create( Params &params, Patch *patch )
    {
    
        std::vector<EnvelopeBC *> EnvBoundCond;
        
        // periodic (=NULL) boundary conditions
        EnvBoundCond.resize( 2*params.nDim_field, NULL );
        
        // -----------------
        // For 2D and 3Dcartesian Geometry only for the moment
        // -----------------
        if( params.geometry == "1Dcartesian" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
                // X DIRECTION
                // reflective bcs
                if( params.Env_BCs[0][ii] == "reflective" ) {
                    EnvBoundCond[ii] = new EnvelopeBC1D_refl( params, patch, ii );
                }
                // else: error
                else {
                    ERROR( "Unknown Envelope x-boundary condition `" << params.Env_BCs[0][ii] << "`" );
                }
                
            }
            
        } else if( params.geometry == "2Dcartesian" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
                // X DIRECTION
                // reflective bcs
                if( params.Env_BCs[0][ii] == "reflective" ) {
                    EnvBoundCond[ii] = new EnvelopeBC2D_refl( params, patch, ii );
                }
                // pml bcs
                else if( params.Env_BCs[0][ii] == "PML" ) {
                    EnvBoundCond[ii] = new EnvelopeBC2D_PML( params, patch, ii );
                }
                // else: error
                else {
                    ERROR( "Unknown Envelope x-boundary condition `" << params.Env_BCs[0][ii] << "`" );
                }
                
                // Y DIRECTION
                // reflective bcs
                if( params.Env_BCs[1][ii] == "reflective" ) {
                    EnvBoundCond[ii+2] = new EnvelopeBC2D_refl( params, patch, ii+2 );
                }
                // pml bcs
                else if( params.Env_BCs[1][ii] == "PML" ) {
                    EnvBoundCond[ii+2] = new EnvelopeBC2D_PML( params, patch, ii+2 );
                }
                // else: error
                else {
                    ERROR( "Unknown Envelope y-boundary condition `" << params.Env_BCs[1][ii] << "`" );
                }
                
            }
            
        }  else if( params.geometry == "3Dcartesian" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
                // X DIRECTION
                // reflective bcs
                if( params.Env_BCs[0][ii] == "reflective" ) {
                    EnvBoundCond[ii] = new EnvelopeBC3D_refl( params, patch, ii );
                }
                // pml bcs
                else if( params.Env_BCs[0][ii] == "PML" ) {
                    EnvBoundCond[ii] = new EnvelopeBC3D_PML( params, patch, ii );
                }
                // else: error
                else {
                    ERROR( "Unknown Envelope x-boundary condition `" << params.Env_BCs[0][ii] << "`" );
                }
                
                // Y DIRECTION
                // reflective bcs
                if( params.Env_BCs[1][ii] == "reflective" ) {
                    EnvBoundCond[ii+2] = new EnvelopeBC3D_refl( params, patch, ii+2 );
                }
                // pml bcs
                else if( params.Env_BCs[1][ii] == "PML" ) {
                    EnvBoundCond[ii+2] = new EnvelopeBC3D_PML( params, patch, ii+2 );
                }
                // else: error
                else {
                    ERROR( "Unknown Envelope y-boundary condition `" << params.Env_BCs[1][ii] << "`" );
                }
                
                // Z DIRECTION
                // reflective bcs
                if( params.Env_BCs[2][ii] == "reflective" ) {
                    EnvBoundCond[ii+4] = new EnvelopeBC3D_refl( params, patch, ii+4 );
                }
                // pml bcs
                else if( params.Env_BCs[2][ii] == "PML" ) {
                    EnvBoundCond[ii+4] = new EnvelopeBC3D_PML( params, patch, ii+4 );
                }
                // else: error
                else  {
                    ERROR( "Unknown Envelope z-boundary condition `" << params.Env_BCs[2][ii] << "`" );
                }
            }
            
        }  else if( params.geometry == "AMcylindrical" ) {

            for( unsigned int ii=0; ii<2; ii++ ) {
            
                // X DIRECTION
                if( params.Env_BCs[0][ii] == "reflective" ) {
                    EnvBoundCond[ii] = new EnvelopeBCAM_refl( params, patch, ii );
                }
                // pml bcs
                else if( params.Env_BCs[0][ii] == "PML" ) {
                    EnvBoundCond[ii] = new EnvelopeBCAM_PML( params, patch, ii );
                } 
                else {
                    ERROR( "Unknown Envelope x-boundary condition `" << params.Env_BCs[0][ii] << "`" );
                }
            }
            
            // R DIRECTION
            EnvBoundCond[2] = new EnvelopeBCAM_Axis( params, patch, 2 );
            
          
            if( params.Env_BCs[1][1] == "reflective" ) {
                EnvBoundCond[3] = new EnvelopeBCAM_refl( params, patch, 3 );  
            }
            // pml bcs
            else if( params.Env_BCs[1][1] == "PML" ) {
                EnvBoundCond[3] = new EnvelopeBCAM_PML( params, patch, 3 );
            }
            else  {
                ERROR( "Unknown Envelope r-boundary condition `" << params.Env_BCs[1][1] << "`" );
            }
        }
        
        
        // OTHER GEOMETRIES ARE NOT DEFINED ---
        else {
            ERROR( "Unknown geometry : " << params.geometry );
        }
        
        return EnvBoundCond;
    }
    
};

#endif

