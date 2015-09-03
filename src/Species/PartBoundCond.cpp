#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include <cmath>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( PicParams& params, int ispec, SmileiMPI* smpi )
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
    bc_west   = NULL;
    bc_east   = NULL;
    bc_south  = NULL;
    bc_north  = NULL;
    bc_bottom = NULL;
    bc_up     = NULL;
    
    // -----------------------------
    // Define limits of local domain
    // -----------------------------
    
    // 1d3v or 2d3v or 3d3v
    x_min = max( x_min_global, smpi->getDomainLocalMin(0) );
    x_max = min( x_max_global, smpi->getDomainLocalMax(0) );
    
    // 2d3v or 3d3v
    if ( nDim_particle > 1 ) {
        if ( (params.bc_em_type_y[0]=="periodic") || (params.bc_em_type_y[1]=="periodic") ) {
            y_min = smpi->getDomainLocalMin(1);
            y_max = smpi->getDomainLocalMax(1);
        }
        else {
            y_min = max( y_min_global, smpi->getDomainLocalMin(1) );
            y_max = min( y_max_global, smpi->getDomainLocalMax(1) );
        }
    }
    
    // 3d3v
    if ( nDim_particle > 2 ) {
        if ( (params.bc_em_type_z[0]=="periodic") || (params.bc_em_type_z[1]=="periodic") ) {
            z_min = smpi->getDomainLocalMin(2);
            z_max = smpi->getDomainLocalMax(2);
        }
        else {
            z_min = max( z_min_global, smpi->getDomainLocalMin(2) );
            z_max = min( z_max_global, smpi->getDomainLocalMax(2) );
        }
    }
    
    // Check for inconsistencies between EM and particle BCs
    if ( ((params.bc_em_type_x[0]=="periodic")&&(params.species_param[ispec].bc_part_type_west!="none"))
     ||  ((params.bc_em_type_x[1]=="periodic")&&(params.species_param[ispec].bc_part_type_east!="none")) ) {
        ERROR("For species #" << ispec << ", periodic EM boundary conditions require x particle BCs to be periodic.");
    }
    if ( nDim_particle > 1 ) {
        if ( ((params.bc_em_type_y[0]=="periodic")&&(params.species_param[ispec].bc_part_type_south!="none"))
         ||  ((params.bc_em_type_y[1]=="periodic")&&(params.species_param[ispec].bc_part_type_north!="none")) ) {
            ERROR("For species #" << ispec << ", periodic EM boundary conditions require y particle BCs to be periodic.");
        }
        if ( nDim_particle > 2 ) {
            if ( ((params.bc_em_type_z[0]=="periodic")&&(params.species_param[ispec].bc_part_type_bottom!="none"))
             ||  ((params.bc_em_type_z[1]=="periodic")&&(params.species_param[ispec].bc_part_type_up!="none"    )) ) {
                ERROR("For species #" << ispec << ", periodic EM boundary conditions require z particle BCs to be periodic.");
            }
        }
    }
    
    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------
    
    bool thermCond = false;
    
    // West
    if ( params.species_param[ispec].bc_part_type_west == "refl" ) {
        if (smpi->isWestern()) bc_west = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "supp" ) {
        if (smpi->isWestern()) bc_west = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "stop" ) {
        if (smpi->isWestern()) bc_west = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "thermalize" ) {
        thermCond = true;
        if (smpi->isWestern()) bc_west = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "none" ) {
        WARNING( "West boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "West boundary condition undefined" );
    }
    
    // East
    if ( params.species_param[ispec].bc_part_type_east == "refl" ) {
        if (smpi->isEastern()) bc_east = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "supp" ) {
        if (smpi->isEastern()) bc_east = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "stop" ) {
        if (smpi->isEastern()) bc_east = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "thermalize" ) {
        thermCond = true;
        if (smpi->isEastern()) bc_east = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "none" ) {
        WARNING( "East boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "East boundary condition undefined" );
    }
    
    
    if ( nDim_particle > 1 ) {
        // South
        if ( params.species_param[ispec].bc_part_type_south == "refl" ) {
            if (smpi->isSouthern()) bc_south = &refl_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "supp" ) {
            if (smpi->isSouthern()) bc_south = &supp_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "stop" ) {
            if (smpi->isSouthern()) bc_south = &stop_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "thermalize" ) {
            thermCond = true;
            if (smpi->isSouthern()) bc_south = &thermalize_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_south == "none" ) {
            WARNING( "South boundary condition for species " << ispec << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "South boundary condition undefined : " << params.species_param[ispec].bc_part_type_south  );
        }
        
        // North
        if ( params.species_param[ispec].bc_part_type_north == "refl" ) {
            if (smpi->isNorthern()) bc_north = &refl_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "supp" ) {
            if (smpi->isNorthern()) bc_north = &supp_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "stop" ) {
            if (smpi->isNorthern()) bc_north = &stop_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "thermalize" ) {
            thermCond = true;
            if (smpi->isNorthern()) bc_north = &thermalize_particle;
        }
        else if ( params.species_param[ispec].bc_part_type_north == "none" ) {
            WARNING( "North boundary condition for species " << ispec << " is 'none', which means the same as fields");
        }
        else {
            ERROR( "North boundary condition undefined : " << params.species_param[ispec].bc_part_type_north  );
        }
        
        
        if ( nDim_particle > 2 ) {
            if ( params.species_param[ispec].bc_part_type_bottom == "refl" ) {
                if (z_min==z_min_global) bc_bottom = &refl_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_bottom == "supp" ) {
                if (z_min==z_min_global) bc_bottom = &supp_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_bottom == "stop" ) {
                if (z_min==z_min_global) bc_bottom = &stop_particle;
            }
            
            if ( params.species_param[ispec].bc_part_type_up == "refl" ) {
                if (z_min==z_min_global) bc_up = &refl_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_up == "supp" )  {
                if (z_min==z_min_global) bc_up = &supp_particle;
            }
            else if ( params.species_param[ispec].bc_part_type_up == "stop" ) {
                if (z_min==z_min_global) bc_up = &stop_particle;
            }
            
        }//nDim_particle>2
        
    }//nDim_particle>1
    
    /* NOT USED ANYMORE AS WE USE THE ERFINV FCT FROM TOOLS/USERFUNCTIONS
    // ---------------------------------------------------------------------
    // Compute the tabulated inverse error function used in thermalizing bcs
    // ---------------------------------------------------------------------
    if ( thermCond ) {
        erfinv::instance().prepare();
    }//thermCond
     */
    
}


PartBoundCond::~PartBoundCond()
{
}


void PartBoundCond::moveWindow_x(double shift, SmileiMPI* smpi)
{
    x_min += shift;
    x_max += shift;
    if (smpi->isWestern()) bc_west = &supp_particle;
    if (smpi->isEastern()) bc_east = &supp_particle;
}
