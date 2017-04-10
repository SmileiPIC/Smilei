
#ifndef ELECTROMAGNBC_FACTORY_H
#define ELECTROMAGNBC_FACTORY_H

#include "ElectroMagnBC.h"
#include "ElectroMagnBC1D_SM.h"
#include "ElectroMagnBC1D_refl.h"
#include "ElectroMagnBC2D_SM.h"
#include "ElectroMagnBC2D_refl.h"
#include "ElectroMagnBC3D_SM.h"

#include "Params.h"


//  --------------------------------------------------------------------------------------------------------------------
//! Constructor for the ElectroMagnetic Boundary conditions factory
//  --------------------------------------------------------------------------------------------------------------------
class ElectroMagnBC_Factory {
    
public:
    
    static std::vector<ElectroMagnBC*> create(Params& params, Patch* patch) {
        
        std::vector<ElectroMagnBC*> emBoundCond;

        // periodic (=NULL) boundary conditions
        emBoundCond.resize(2*params.nDim_field, NULL);
        
        // -----------------
        // For 1d3v Geometry
        // -----------------
        if ( params.geometry == "1d3v" ) {
            
            
            // AT X = XMIN,XMAX
            // ----------------
            for (unsigned int ii=0;ii<2;ii++) {
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_SM(params, patch, ii);
                }
                // reflective bcs
                else if ( params.bc_em_type_x[ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_refl(params, patch, ii);
                }
                // else: error
                else if ( params.bc_em_type_x[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
                }
            }
            
        }//1d3v
        
        
        // -----------------
        // For 2d3v Geometry
        // -----------------
        else if ( params.geometry == "2d3v" ) {
            
            for (unsigned int ii=0;ii<2;ii++) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_SM(params, patch, ii);
                }
                // reflective bcs
                else if ( params.bc_em_type_x[ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_refl(params, patch, ii);
                }
                // else: error
                else if ( params.bc_em_type_x[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_y[ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_SM(params, patch, ii+2);
                }
                // reflective bcs
                else if ( params.bc_em_type_y[ii] == "reflective" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_refl(params, patch, ii+2);
                }
                // else: error
                else if ( params.bc_em_type_y[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }
            }
            
        }//2d3v
        
        // -----------------
        // For 3d3v Geometry
        // -----------------
        else if ( params.geometry == "3d3v" ) {
            
            for (unsigned int ii=0;ii<2;ii++) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC3D_SM(params, patch, ii);
                }
                // else: error
                else if ( params.bc_em_type_x[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_y[ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC3D_SM(params, patch, ii+2);
                }
                // else: error
                else if ( params.bc_em_type_y[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }

                // Z DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_z[ii] == "silver-muller" ) {
                    emBoundCond[ii+4] = new ElectroMagnBC3D_SM(params, patch, ii+4);
                }
                // else: error
                else if ( params.bc_em_type_z[ii] != "periodic" ) {
                    ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }
            }
            
        }//3d3v       


        // OTHER GEOMETRIES ARE NOT DEFINED ---
        else {
            ERROR( "Unknown geometry : " << params.geometry );
        }
        
        return emBoundCond;
    }
    
};

#endif

