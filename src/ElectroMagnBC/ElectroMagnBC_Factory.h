
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
        
        unsigned int nperiodicx=0, nperiodicy=0, nperiodicz=0;
        
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
                    emBoundCond[ii] = new ElectroMagnBC1D_SM(params, patch);
                }
                // reflective bcs
                else if ( params.bc_em_type_x[ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC1D_refl(params, patch);
                }
                // else: error
                else {
                    nperiodicx++;
                    if ( params.bc_em_type_x[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
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
                    emBoundCond[ii] = new ElectroMagnBC2D_SM(params, patch);
                }
                // reflective bcs
                else if ( params.bc_em_type_x[ii] == "reflective" ) {
                    emBoundCond[ii] = new ElectroMagnBC2D_refl(params, patch);
                }
                // else: error
                else {
                    nperiodicx++;
                    if ( params.bc_em_type_x[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_y[ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_SM(params, patch);
                }
                // reflective bcs
                else if ( params.bc_em_type_y[ii] == "reflective" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC2D_refl(params, patch);
                }
                // else: error
                else {
                    nperiodicy++;
                    if ( params.bc_em_type_y[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }
            }
            
        }//2d3v
        
        // -----------------
        // For 3d3v Geometry
        // -----------------
        else if ( params.geometry == "3d3v" ) {
            
            // by default use periodic (=NULL) boundary conditions
            emBoundCond.resize(6, NULL);
            
            for (unsigned int ii=0;ii<2;ii++) {
                // X DIRECTION
                // silver-muller (injecting/absorbing bcs)
                if ( params.bc_em_type_x[ii] == "silver-muller" ) {
                    emBoundCond[ii] = new ElectroMagnBC3D_SM(params, patch);
                }
                // else: error
                else {
                    nperiodicx++;
                    if ( params.bc_em_type_x[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_x[" << ii << "]");
                }
                
                // Y DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_y[ii] == "silver-muller" ) {
                    emBoundCond[ii+2] = new ElectroMagnBC3D_SM(params, patch);
                }
                // else: error
                else {
                    nperiodicy++;
                    if ( params.bc_em_type_y[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }
                
                // Z DIRECTION
                // silver-muller bcs (injecting/absorbin)
                if ( params.bc_em_type_z[ii] == "silver-muller" ) {
                    emBoundCond[ii+4] = new ElectroMagnBC3D_SM(params, patch);
                }
                // else: error
                else {
                    nperiodicz++;
                    if ( params.bc_em_type_z[ii] != "periodic" )
                        ERROR( "Unknown boundary bc_em_type_y[" << ii << "]");
                }
            }
            
        }//3d3v
        
        
        // OTHER GEOMETRIES ARE NOT DEFINED ---
        else {
            ERROR( "Unknown geometry : " << params.geometry );
        }
        
        
        if( nperiodicx == 1 )
            ERROR("EM boundary conditions, in direction x, cannot be periodic on one side and not on the other");
        if( nperiodicy == 1 )
            ERROR("EM boundary conditions, in direction y, cannot be periodic on one side and not on the other");
        if( nperiodicz == 1 )
            ERROR("EM boundary conditions, in direction z, cannot be periodic on one side and not on the other");
        
        return emBoundCond;
    }
    
};

#endif

