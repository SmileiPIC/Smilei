
#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC1D_refl.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC2D_refl.h"
#include "ElectroMagnBC3D_SM.h"
#include "ElectroMagnBC3D_refl.h"
#include "ElectroMagnBC3D_BM.h"
#include "ElectroMagnBCAM_SM.h"
#include "ElectroMagnBCAM_zero.h"
#include "ElectroMagnBCAM_BM.h"

#include "Params.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Constructor for the ElectroMagnetic Boundary conditions factory
//  --------------------------------------------------------------------------------------------------------------------
class ElectroMagnBC_Factory
{

public:

    static std::vector<ElectroMagnBC *> create( Params &params, Patch *patch )
    {
    
        std::vector<ElectroMagnBC *> emBoundCond;
        
        // periodic (=NULL) boundary conditions
        emBoundCond.resize( 2*params.nDim_field, NULL );
        
        // -----------------
        // For 1Dcartesian Geometry
        // -----------------
        if( params.geometry == "1Dcartesian" ) {
        
        
            // AT X = XMIN,XMAX
            // ----------------
            for( unsigned int ii=0; ii<2; ii++ ) {
                // silver-muller (injecting/absorbing bcs)
                if( params.EM_BCs[0][ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_SM( params, patch, ii );
                }
                // reflective bcs
                else if( params.EM_BCs[0][ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_refl( params, patch, ii );
                }
                // else: error
                else if( params.EM_BCs[0][ii] != "periodic" ) {
                    ERROR( "Unknown EM x-boundary condition `" << params.EM_BCs[0][ii] << "`" );
                }
            }
            
        }//1Dcartesian
        
        
        // -----------------
        // For 2Dcartesian Geometry
        // -----------------
        else if( params.geometry == "2Dcartesian" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if( params.EM_BCs[0][ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_SM( params, patch, ii );
                }
                // reflective bcs
                else if( params.EM_BCs[0][ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_refl( params, patch, ii );
                }
                // else: error
                else if( params.EM_BCs[0][ii] != "periodic" ) {
                    ERROR( "Unknown EM x-boundary condition `" << params.EM_BCs[0][ii] << "`" );
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if( params.EM_BCs[1][ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_SM( params, patch, ii+2 );
                }
                // reflective bcs
                else if( params.EM_BCs[1][ii] == "reflective" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_refl( params, patch, ii+2 );
                }
                // else: error
                else if( params.EM_BCs[1][ii] != "periodic" ) {
                    ERROR( "Unknown EM y-boundary condition `" << params.EM_BCs[1][ii] << "`" );
                }
            }
            
        }//2Dcartesian
        
        // -----------------
        // For 3Dcartesian Geometry
        // -----------------
        else if( params.geometry == "3Dcartesian" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if( params.EM_BCs[0][ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC3D_SM( params, patch, ii );
                }
                // reflective bcs
                else if( params.EM_BCs[0][ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC3D_refl( params, patch, ii );
                }
                // Buneman bcs (absorbing)
                else if( params.EM_BCs[0][ii] == "buneman" ) {
                    emBoundCond[ii] = new ElectroMagnBC3D_BM( params, patch, ii );
                }
                // else: error
                else if( params.EM_BCs[0][ii] != "periodic" ) {
                    ERROR( "Unknown EM x-boundary condition `" << params.EM_BCs[0][ii] << "`" );
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbing)
                if( params.EM_BCs[1][ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC3D_SM( params, patch, ii+2 );
                }
                // reflective bcs
                else if( params.EM_BCs[1][ii] == "reflective" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC3D_refl( params, patch, ii+2 );
                }
                // Buneman bcs (absorbing)
                else if( params.EM_BCs[1][ii] == "buneman" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC3D_BM( params, patch, ii+2 );
                }
                // else: error
                else if( params.EM_BCs[1][ii] != "periodic" ) {
                    ERROR( "Unknown EM y-boundary condition `" << params.EM_BCs[1][ii] << "`" );
                }
                
                // Z DIRECTION
                // silver-muller bcs (injecting/absorbing)
                if( params.EM_BCs[2][ii] == "silver-muller" ) {
                    emBoundCond[ii+4] = new ElectroMagnBC3D_SM( params, patch, ii+4 );
                }
                // reflective bcs
                else if( params.EM_BCs[2][ii] == "reflective" ) {
                    emBoundCond[ii+4] = new ElectroMagnBC3D_refl( params, patch, ii+4 );
                }
                // Buneman bcs (absorbing)
                else if( params.EM_BCs[2][ii] == "buneman" ) {
                    emBoundCond[ii+4] = new ElectroMagnBC3D_BM( params, patch, ii+4 );
                }
                // else: error
                else if( params.EM_BCs[2][ii] != "periodic" ) {
                    ERROR( "Unknown EM z-boundary condition `" << params.EM_BCs[2][ii] << "`" );
                }
            }
            
        }//3Dcartesian
        
        // -----------------
        // For AM cylindrical Geometry
        // -----------------
        else if( params.geometry == "AMcylindrical" ) {
        
            for( unsigned int ii=0; ii<2; ii++ ) {
            
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if( params.EM_BCs[0][ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBCAM_SM( params, patch, ii );
                }
                else if( params.EM_BCs[0][ii] == "zero" ) {
                    emBoundCond[ii] = new ElectroMagnBCAM_zero( params, patch, ii );
                }
                else if( params.EM_BCs[0][ii] != "periodic" ) {
                    ERROR( "Unknown EM x-boundary condition `" << params.EM_BCs[0][ii] << "`" );
                }
            }
            
            // R DIRECTION
            emBoundCond[2] = NULL ; //Axis BC are handeled directly in solvers.
            if( params.EM_BCs[1][1] == "periodic" ) {
                ERROR( "Periodic EM Rmax-boundary condition is not supported`" );
            } else if( params.EM_BCs[1][1] == "buneman" ) {
                emBoundCond[3] = new ElectroMagnBCAM_BM( params, patch, 3 );
            } else if (params.is_spectral) {
                emBoundCond[3] = NULL ;
            } else  {
                ERROR( "Unknown EM Rmax-boundary condition `" << params.EM_BCs[1][1] << "`" );
            }
            //MESSAGE( params.EM_BCs[1][1]);
            
        }//AM
        
        // OTHER GEOMETRIES ARE NOT DEFINED ---
        else {
            ERROR( "Unknown geometry : " << params.geometry );
        }
        
        return emBoundCond;
    }
    
};

#endif

