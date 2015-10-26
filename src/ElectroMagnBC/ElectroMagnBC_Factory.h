
#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC1D_refl.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC2D_refl.h"

#include "Params.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Constructor for the ElectroMagnetic Boundary conditions factory
//  --------------------------------------------------------------------------------------------------------------------
class ElectroMagnBC_Factory {
    
public:
    
    static std::vector<ElectroMagnBC*> create(Params& params, LaserParams &laser_params) {
        
        std::vector<ElectroMagnBC*> emBoundCond;
        
        // -----------------
        // For 1d3v Geometry
        // -----------------
        if ( params.geometry == "1d3v" ) {
            
            // periodic (=NULL) boundary conditions
            emBoundCond.resize(2, NULL);
            
            // AT X = XMIN,XMAX
            // ----------------
            for (unsigned int ii=0;ii<2;ii++) {
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_SM(params, laser_params);
                }
                // reflective bcs
                else if ( params.bc_em_type_x[i] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_refl(params, laser_params);
                }
                // else: error
                else if ( params.bc_em_type_x[ii] != "periodic" ) {
                    ERROR( "Unknwon boundary bc_em_type_x[" << ii << "]");
                }
            }
            
        }//1d3v
        
        
        // -----------------
        // For 2d3v Geometry
        // -----------------
        else if ( params.geometry == "2d3v" ) {
            
            // by default use periodic (=NULL) boundary conditions
            emBoundCond.resize(4, NULL);

            for (unsigned int ii=0;ii<2;ii++) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_SM(params, laser_params);
                }
                // reflective bcs
                else if ( params.   [i] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_refl(params, laser_params);
                }
                // else: error
                else if ( params.bc_em_type_x[ii] != "periodic" ) {
                    ERROR( "Unknwon boundary bc_em_type_x[" << ii << "]");
                }

                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_y[0] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_SM(params, laser_params);
                }
                // reflective bcs
                else if ( params.bc_em_type_y[0] == "reflective" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_refl(params, laser_params);
                }
                // else: error
                else if ( params.bc_em_type_y[0] != "periodic" ) {
                    ERROR( "Unknwon boundary bc_em_type_y[" << ii << "]");
                }
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

