
#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC1D_refl.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC2D_refl.h"

#include "PicParams.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Constructor for the ElectroMagnetic Boundary conditions factory
//  --------------------------------------------------------------------------------------------------------------------
class ElectroMagnBC_Factory {
    
public:
    
    static std::vector<ElectroMagnBC*> create(PicParams& params, LaserParams &laser_params) {
        
        std::vector<ElectroMagnBC*> emBoundCond;
        
        // -----------------
        // For 1d3v Geometry
        // -----------------
        if ( params.geometry == "1d3v" ) {
            
            // periodic (=NULL) boundary conditions
            emBoundCond.resize(2, NULL);
            
            // AT X = XMIN
            // -----------
            
            // silver-muller (injecting/absorbing bcs)
            if ( params.bc_em_type_x[0] == "silver-muller" ) {
                emBoundCond[0] = new ElectroMagnBC1D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_x[0] == "reflective" ) {
                emBoundCond[0] = new ElectroMagnBC1D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_x[0] != "periodic" ) {
                ERROR( "Unknwon boundary condition at x=xmin: " << params.bc_em_type_x[0] );
            }
            
            // AT X = XMAX
            // -----------
            
            // silver-muller (injecting/absorbing bcs)
            if ( params.bc_em_type_x[1] == "silver-muller" ) {
                emBoundCond[1] = new ElectroMagnBC1D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_x[1] == "reflective" ) {
                //emBoundCond[1] = new ElectroMagnBC1D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_x[1] != "periodic" ) {
                ERROR( "Unknwon boundary condition at x=xmax: " << params.bc_em_type_x[1] );
            }
            
        }//1d3v
        
        
        // -----------------
        // For 2d3v Geometry
        // -----------------
        else if ( params.geometry == "2d3v" ) {
            
            // by default use periodic (=NULL) boundary conditions
            emBoundCond.resize(4, NULL);

            // X-Direction: @ xmin
            // -------------------
            
            // silver-muller bcs (injecting/absorbin)
            if ( params.bc_em_type_x[0] == "silver-muller" ) {
                emBoundCond[0] = new ElectroMagnBC2D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_x[0] == "reflective" ) {
                emBoundCond[0] = new ElectroMagnBC2D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_x[0] != "periodic" ) {
                ERROR( "Unknwon boundary condition at x=xmin: " << params.bc_em_type_x[0] );
            }
            
            // X-Direction: @ xmax
            // -------------------
            
            // silver-muller bcs (injecting/absorbin)
            if ( params.bc_em_type_x[1] == "silver-muller" ) {
                emBoundCond[1] = new ElectroMagnBC2D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_x[1] == "reflective" ) {
                emBoundCond[1] = new ElectroMagnBC2D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_x[1] != "periodic" ) {
                ERROR( "Unknwon boundary condition at x=xmax: " << params.bc_em_type_x[0] );
            }
            
            // Y-Direction: @ ymin
            // -------------------
            
            // silver-muller bcs (injecting/absorbin)
            if ( params.bc_em_type_y[0] == "silver-muller" ) {
                emBoundCond[2] = new ElectroMagnBC2D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_y[0] == "reflective" ) {
                emBoundCond[2] = new ElectroMagnBC2D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_y[0] != "periodic" ) {
                ERROR( "Unknwon boundary condition at y=ymin: " << params.bc_em_type_y[0] );
            }
            
            // Y-Direction: @ ymax
            // -------------------
            
            // silver-muller bcs (injecting/absorbin)
            if ( params.bc_em_type_y[1] == "silver-muller" ) {
                emBoundCond[3] = new ElectroMagnBC2D_SM(params, laser_params);
            }
            // reflective bcs
            else if ( params.bc_em_type_y[1] == "reflective" ) {
                emBoundCond[3] = new ElectroMagnBC2D_refl(params, laser_params);
            }
            // else: error
            else if ( params.bc_em_type_y[1] != "periodic" ) {
                ERROR( "Unknwon boundary condition at y=ymax: " << params.bc_em_type_y[1] );
            }

        }//2d3v
        
        
        // OTHER GEOMETRIES ARE NOT DEFINED ---
        else {
            ERROR( "Unknwon geometry : " << params.geometry );
        }
        return emBoundCond;
    }
    
};

#endif

