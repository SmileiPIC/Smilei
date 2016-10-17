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
    if (params.bc_em_type_x[0]=="periodic" || params.hasWindow) {
        x_min = patch->getDomainLocalMin(0);
        x_max = patch->getDomainLocalMax(0);
    }
    else {
        x_min = max( x_min_global, patch->getDomainLocalMin(0) );
        x_max = min( x_max_global, patch->getDomainLocalMax(0) );
    }
    
    if ( nDim_particle > 1 ) {
        if (params.bc_em_type_y[0]=="periodic") {
            y_min = patch->getDomainLocalMin(1);
            y_max = patch->getDomainLocalMax(1);
        }
        else {
            y_min = max( y_min_global, patch->getDomainLocalMin(1) );
            y_max = min( y_max_global, patch->getDomainLocalMax(1) );
        }
        
        if ( nDim_particle > 2 ) {
            if (params.bc_em_type_z[0]=="periodic") {
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
        if ( ((params.bc_em_type_x[0]=="periodic")&&(species->bc_part_type_xmin!="none"))
         ||  ((params.bc_em_type_x[1]=="periodic")&&(species->bc_part_type_xmax!="none")) ) {
            ERROR("For species " << species->species_type << ", periodic EM boundary conditions require x particle BCs to be periodic.");
        }
        if ( nDim_particle > 1 ) {
            if ( ((params.bc_em_type_y[0]=="periodic")&&(species->bc_part_type_ymin!="none"))
             ||  ((params.bc_em_type_y[1]=="periodic")&&(species->bc_part_type_ymax!="none")) ) {
                ERROR("For species #" << species->species_type << ", periodic EM boundary conditions require y particle BCs to be periodic.");
            }
            if ( nDim_particle > 2 ) {
                if ( ((params.bc_em_type_z[0]=="periodic")&&(species->bc_part_type_zmin!="none"))
                 ||  ((params.bc_em_type_z[1]=="periodic")&&(species->bc_part_type_zmax!="none"    )) ) {
                    ERROR("For species #" << species->species_type << ", periodic EM boundary conditions require z particle BCs to be periodic.");
                }
            }
        }
    }
    
    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------
    
    // West
    if ( species->bc_part_type_xmin == "refl" ) {
        if (patch->isXmin()) bc_xmin = &refl_particle;
    }
    else if ( species->bc_part_type_xmin == "supp" ) {
        if (patch->isXmin()) bc_xmin = &supp_particle;
    }
    else if ( species->bc_part_type_xmin == "stop" ) {
        if (patch->isXmin()) bc_xmin = &stop_particle;
    }
    else if ( species->bc_part_type_xmin == "thermalize" ) {
        if (patch->isXmin()) bc_xmin = &thermalize_particle;
    }
    else if ( species->bc_part_type_xmin == "none" ) {
        if (patch->isMaster()) MESSAGE(2,"West boundary condition for species " << species->species_type << " is 'none', which means the same as fields");
    }
    else {
        ERROR("West boundary condition undefined" );
    }
    
    // East
    if ( species->bc_part_type_xmax == "refl" ) {
        if (patch->isXmax()) bc_xmax = &refl_particle;
    }
    else if ( species->bc_part_type_xmax == "supp" ) {
        if (patch->isXmax()) bc_xmax = &supp_particle;
    }
    else if ( species->bc_part_type_xmax == "stop" ) {
        if (patch->isXmax()) bc_xmax = &stop_particle;
    }
    else if ( species->bc_part_type_xmax == "thermalize" ) {
        if (patch->isXmax()) bc_xmax = &thermalize_particle;
    }
    else if ( species->bc_part_type_xmax == "none" ) {
        if (patch->isMaster()) MESSAGE(2,"East boundary condition for species " << species->species_type << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "East boundary condition undefined" );
    }
    
    
    if ( nDim_particle > 1 ) {
        // South
        if ( species->bc_part_type_ymin == "refl" ) {
            if (patch->isYmin()) bc_ymin = &refl_particle;
        }
        else if ( species->bc_part_type_ymin == "supp" ) {
            if (patch->isYmin()) bc_ymin = &supp_particle;
        }
        else if ( species->bc_part_type_ymin == "stop" ) {
            if (patch->isYmin()) bc_ymin = &stop_particle;
        }
        else if ( species->bc_part_type_ymin == "thermalize" ) {
            if (patch->isYmin()) bc_ymin = &thermalize_particle;
        }
        else if ( species->bc_part_type_ymin == "none" ) {
            if (patch->isMaster()) MESSAGE(2,"South boundary condition for species " << species->species_type << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "South boundary condition undefined : " << species->bc_part_type_ymin  );
        }
        
        // North
        if ( species->bc_part_type_ymax == "refl" ) {
            if (patch->isYmax()) bc_ymax = &refl_particle;
        }
        else if ( species->bc_part_type_ymax == "supp" ) {
            if (patch->isYmax()) bc_ymax = &supp_particle;
        }
        else if ( species->bc_part_type_ymax == "stop" ) {
            if (patch->isYmax()) bc_ymax = &stop_particle;
        }
        else if ( species->bc_part_type_ymax == "thermalize" ) {
            if (patch->isYmax()) bc_ymax = &thermalize_particle;
        }
        else if ( species->bc_part_type_ymax == "none" ) {
            if (patch->isMaster()) MESSAGE(2,"North boundary condition for species " << species->species_type << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "North boundary condition undefined : " << species->bc_part_type_ymax  );
        }
        
        
        if ( nDim_particle > 2 ) {
            if ( species->bc_part_type_zmin == "refl" ) {
                if (z_min==z_min_global) bc_zmin = &refl_particle;
            }
            else if ( species->bc_part_type_zmin == "supp" ) {
                if (z_min==z_min_global) bc_zmin = &supp_particle;
            }
            else if ( species->bc_part_type_zmin == "stop" ) {
                if (z_min==z_min_global) bc_zmin = &stop_particle;
            }
            
            if ( species->bc_part_type_zmax == "refl" ) {
                if (z_min==z_min_global) bc_zmax = &refl_particle;
            }
            else if ( species->bc_part_type_zmax == "supp" )  {
                if (z_min==z_min_global) bc_zmax = &supp_particle;
            }
            else if ( species->bc_part_type_zmax == "stop" ) {
                if (z_min==z_min_global) bc_zmax = &stop_particle;
            }
            
        }//nDim_particle>2
        
    }//nDim_particle>1
    
    
}


PartBoundCond::~PartBoundCond()
{
}

