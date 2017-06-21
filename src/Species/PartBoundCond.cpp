#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( Params& params, Species *species, Patch* patch )
{
    // number of dimensions for the particle
    //!\todo (MG to JD) isn't it always 3?
    nDim_particle = params.nDim_particle;
    
    // Absolute global values
    double x_min_global = 0;
    double x_max_global = params.cell_length[0]*(params.n_space_global[0]);
    double y_min_global = 0;
    double y_max_global = params.cell_length[1]*(params.n_space_global[1]);
    double z_min_global = 0;
    double z_max_global = params.cell_length[2]*(params.n_space_global[2]);
    
    // by default apply no bcs
    bc_xmin   = NULL;
    bc_xmax   = NULL;
    bc_ymin  = NULL;
    bc_ymax  = NULL;
    bc_zmin = NULL;
    bc_zmax     = NULL;
    
    // -----------------------------
    // Define limits of local domain
    if (params.EM_BCs[0][0]=="periodic" || params.hasWindow) {
        x_min = patch->getDomainLocalMin(0);
        x_max = patch->getDomainLocalMax(0);
    }
    else {
        x_min = max( x_min_global, patch->getDomainLocalMin(0) );
        x_max = min( x_max_global, patch->getDomainLocalMax(0) );
    }
    
    if ( nDim_particle > 1 ) {
        if (params.EM_BCs[1][0]=="periodic") {
            y_min = patch->getDomainLocalMin(1);
            y_max = patch->getDomainLocalMax(1);
        }
        else {
            y_min = max( y_min_global, patch->getDomainLocalMin(1) );
            y_max = min( y_max_global, patch->getDomainLocalMax(1) );
        }
        
        if ( nDim_particle > 2 ) {
            if (params.EM_BCs[2][0]=="periodic") {
                z_min = patch->getDomainLocalMin(2);
                z_max = patch->getDomainLocalMax(2);
            }
            else {
                z_min = max( z_min_global, patch->getDomainLocalMin(2) );
                z_max = min( z_max_global, patch->getDomainLocalMax(2) );
            }
        }
    }


    // Can be done after parsing 

    // Check for inconsistencies between EM and particle BCs
    if (! species->particles->tracked) {
        for( unsigned int iDim=0; iDim<nDim_particle; iDim++ ) {
            if ( ((params.EM_BCs[iDim][0]=="periodic")&&(species->boundary_conditions[iDim][0]!="periodic"))  
             ||  ((params.EM_BCs[iDim][1]=="periodic")&&(species->boundary_conditions[iDim][1]!="periodic")) ) {
                ERROR("For species " << species->species_type << ", periodic EM "<<"xyz"[iDim]<<"-boundary conditions require particle BCs to be periodic.");
            }
        }
    }
    
    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------
    
    // Xmin
    if ( species->boundary_conditions[0][0] == "refl" ) {
        if (patch->isXmin()) bc_xmin = &refl_particle;
    }
    else if ( species->boundary_conditions[0][0] == "supp" ) {
        if (patch->isXmin()) bc_xmin = &supp_particle;
    }
    else if ( species->boundary_conditions[0][0] == "stop" ) {
        if (patch->isXmin()) bc_xmin = &stop_particle;
    }
    else if ( species->boundary_conditions[0][0] == "thermalize" ) {
        if (patch->isXmin()) bc_xmin = &thermalize_particle;
    }
    else if ( species->boundary_conditions[0][0] == "periodic" ) {
    }
    else {
        ERROR("Xmin boundary condition `"<<species->boundary_conditions[0][0]<<"` unknown" );
    }
    
    // Xmax
    if ( species->boundary_conditions[0][1] == "refl" ) {
        if (patch->isXmax()) bc_xmax = &refl_particle;
    }
    else if ( species->boundary_conditions[0][1] == "supp" ) {
        if (patch->isXmax()) bc_xmax = &supp_particle;
    }
    else if ( species->boundary_conditions[0][1] == "stop" ) {
        if (patch->isXmax()) bc_xmax = &stop_particle;
    }
    else if ( species->boundary_conditions[0][1] == "thermalize" ) {
        if (patch->isXmax()) bc_xmax = &thermalize_particle;
    }
    else if ( species->boundary_conditions[0][1] == "periodic" ) {
    }
    else {
        ERROR( "Xmax boundary condition `"<<species->boundary_conditions[0][1]<<"`  unknown" );
    }
    
    
    if ( nDim_particle > 1 ) {
        // Ymin
        if ( species->boundary_conditions[1][0] == "refl" ) {
            if (patch->isYmin()) bc_ymin = &refl_particle;
        }
        else if ( species->boundary_conditions[1][0] == "supp" ) {
            if (patch->isYmin()) bc_ymin = &supp_particle;
        }
        else if ( species->boundary_conditions[1][0] == "stop" ) {
            if (patch->isYmin()) bc_ymin = &stop_particle;
        }
        else if ( species->boundary_conditions[1][0] == "thermalize" ) {
            if (patch->isYmin()) bc_ymin = &thermalize_particle;
        }
        else if ( species->boundary_conditions[1][0] == "periodic" ) {
        }
        else {
            ERROR( "Ymin boundary condition `"<< species->boundary_conditions[1][0] << "` unknown"  );
        }
        
        // Ymax
        if ( species->boundary_conditions[1][1] == "refl" ) {
            if (patch->isYmax()) bc_ymax = &refl_particle;
        }
        else if ( species->boundary_conditions[1][1] == "supp" ) {
            if (patch->isYmax()) bc_ymax = &supp_particle;
        }
        else if ( species->boundary_conditions[1][1] == "stop" ) {
            if (patch->isYmax()) bc_ymax = &stop_particle;
        }
        else if ( species->boundary_conditions[1][1] == "thermalize" ) {
            if (patch->isYmax()) bc_ymax = &thermalize_particle;
        }
        else if ( species->boundary_conditions[1][1] == "periodic" ) {
        }
        else {
            ERROR( "Ymax boundary condition `"<< species->boundary_conditions[1][1] <<"` undefined" );
        }
        
        
        if ( nDim_particle > 2 ) {
            if ( species->boundary_conditions[2][0] == "refl" ) {
                if (z_min==z_min_global) bc_zmin = &refl_particle;
            }
            else if ( species->boundary_conditions[2][0] == "supp" ) {
                if (z_min==z_min_global) bc_zmin = &supp_particle;
            }
            else if ( species->boundary_conditions[2][0] == "stop" ) {
                if (z_min==z_min_global) bc_zmin = &stop_particle;
            }
            
            if ( species->boundary_conditions[2][1] == "refl" ) {
                if (z_min==z_min_global) bc_zmax = &refl_particle;
            }
            else if ( species->boundary_conditions[2][1] == "supp" )  {
                if (z_min==z_min_global) bc_zmax = &supp_particle;
            }
            else if ( species->boundary_conditions[2][1] == "stop" ) {
                if (z_min==z_min_global) bc_zmax = &stop_particle;
            }
            
        }//nDim_particle>2
        
    }//nDim_particle>1
    
    
}


PartBoundCond::~PartBoundCond()
{
}

