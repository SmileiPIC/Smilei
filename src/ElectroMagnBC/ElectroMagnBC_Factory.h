
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
        HEREIAM(params.geometry);
        std::vector<ElectroMagnBC*> emBoundCond;
        
        // -----------------
        // x direction: all geometries
        // -----------------
        for (unsigned int i=0;i<2;i++) {
            HEREIAM(params.bc_em_type_x[i]);
            // reflective bcs
            if ( params.bc_em_type_x[i] == "periodic" ) {
                emBoundCond.push_back(NULL);
            }
            // reflective bcs
            else if ( params.bc_em_type_x[i] == "silver-muller" ) {
                emBoundCond.push_back(new ElectroMagnBC1D_SM(params, laser_params));
            }
            // reflective bcs
            else if ( params.bc_em_type_x[i] == "reflective" ) {
                emBoundCond.push_back(new ElectroMagnBC1D_refl(params, laser_params));
            }
            // else: error
            else {
                ERROR( "Unknwon boundary condition at x[" << i << "]" << params.bc_em_type_x[i] );
            }
        }
    
        // -----------------
        // y direction: 2d3v and 3d3v
        // -----------------
        if ( params.geometry == "2d3v" || params.geometry == "3d3v") {
            HEREIAM("");
            for (unsigned int i=0;i<2;i++) {
                HEREIAM(params.bc_em_type_y[i]);
                if ( params.bc_em_type_x[i] == "periodic" ) {
                    emBoundCond.push_back(NULL);
                }
                else if ( params.bc_em_type_y[i] == "silver-muller" ) {
                    emBoundCond.push_back(new ElectroMagnBC2D_SM(params, laser_params));
                }
                // reflective bcs
                else if ( params.bc_em_type_y[i] == "reflective" ) {
                    emBoundCond.push_back(new ElectroMagnBC2D_refl(params, laser_params));
                }
                // else: error
                else {
                    ERROR( "Unknwon boundary condition at y[" << i << "]" << params.bc_em_type_y[i] );
                }
            }
        }//

        // -----------------
        // z direction: 3d3v
        // -----------------
        if ( params.geometry == "3d3v") {
            //Write here the 3d part...
        }
        
        HEREIAM(emBoundCond.size());

        return emBoundCond;
    }
    
};

#endif

