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

PartBoundCond::PartBoundCond( Params& params, int ispec, Patch* patch )
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
    if (!params.nspace_win_x) {
        x_min = max( x_min_global, patch->getDomainLocalMin(0) );
        x_max = min( x_max_global, patch->getDomainLocalMax(0) );
    }
    else {
        x_min = patch->getDomainLocalMin(0);
        x_max = patch->getDomainLocalMax(0);
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

    
    // Check for inconsistencies between EM and particle BCs
    if (! params.species_param[ispec].isTest) {
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
    }
    
    // ----------------------------------------------
    // Define the kind of applied boundary conditions
    // ----------------------------------------------
    
    // West
    if ( params.species_param[ispec].bc_part_type_west == "refl" ) {
	if (patch->isWestern()) bc_west = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "supp" ) {
	if (patch->isWestern()) bc_west = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "stop" ) {
	if (patch->isWestern()) bc_west = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "thermalize" ) {
	if (patch->isWestern()) bc_west = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "none" ) {
        MESSAGE( "West boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "West boundary condition undefined" );
    }
    
    // East
    if ( params.species_param[ispec].bc_part_type_east == "refl" ) {
	if (patch->isEastern()) bc_east = &refl_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "supp" ) {
	if (patch->isEastern()) bc_east = &supp_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "stop" ) {
	if (patch->isEastern()) bc_east = &stop_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "thermalize" ) {
	if (patch->isEastern()) bc_east = &thermalize_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "none" ) {
        MESSAGE( "East boundary condition for species " << ispec << " is 'none', which means the same as fields");
    }
    else {
        ERROR( "East boundary condition undefined" );
    }
    
    
    if ( nDim_particle > 1 ) {
	// South 
	if ( params.species_param[ispec].bc_part_type_south == "refl" ) {
	    if (patch->isSouthern()) bc_south = &refl_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_south == "supp" ) {
	    if (patch->isSouthern()) bc_south = &supp_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_south == "stop" ) {
	    if (patch->isSouthern()) bc_south = &stop_particle;
	}	
	else if ( params.species_param[ispec].bc_part_type_south == "thermalize" ) {
	    if (patch->isSouthern()) bc_south = &thermalize_particle;
	}	
	else if ( params.species_param[ispec].bc_part_type_south == "none" ) {
            MESSAGE( "South boundary condition for species " << ispec << " is 'none', which means the same as fields");
	}
	else {
	    ERROR( "South boundary condition undefined : " << params.species_param[ispec].bc_part_type_south  );
	}

	// North
	if ( params.species_param[ispec].bc_part_type_north == "refl" ) {
	    if (patch->isNorthern()) bc_north = &refl_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "supp" ) {
	    if (patch->isNorthern()) bc_north = &supp_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "stop" ) {
	    if (patch->isNorthern()) bc_north = &stop_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "thermalize" ) {
	    if (patch->isNorthern()) bc_north = &thermalize_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "none" ) {
            MESSAGE( "North boundary condition for species " << ispec << " is 'none', which means the same as fields");
	}
	else {
	    ERROR( "North boundary condition undefined : " << params.species_param[ispec].bc_part_type_north  );
	}


	//} // else NULL
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

	    //} // else NULL
            
        }//nDim_particle>2
        
    }//nDim_particle>1
    
    
}


PartBoundCond::~PartBoundCond()
{
}

