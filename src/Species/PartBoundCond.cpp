#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( PicParams& params, int ispec, Patch* patch )
{
    nDim_particle = params.nDim_particle;

    int n_ord_proj_max = 5;
    //!\todo define n_ord_proj_max, value from SQUASH
    if (params.interpolation_order==2) n_ord_proj_max = 5;
    if (params.interpolation_order==3) n_ord_proj_max = 5;

    // Absolute global values
//    double x_min_global = params.cell_length[0]*n_ord_proj_max;
//    double x_max_global = params.cell_length[0]*( params.n_space_global[0]-1-n_ord_proj_max );
    double x_min_global = 0;
    double x_max_global = params.cell_length[0]*(params.n_space_global[0]);

    double y_min_global = 0;
    double y_max_global = params.cell_length[1]*(params.n_space_global[1]);
    double z_min_global = 0;
    double z_max_global = params.cell_length[2]*(params.n_space_global[2]);

    bc_west  =  &supp_particle ;
    bc_east  =  &supp_particle ;
    bc_south =  &supp_particle ;
    bc_north =  &supp_particle ;
    bc_bottom = &supp_particle ;
    bc_up     = &supp_particle ;

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
	if (params.bc_em_type_trans=="periodic") {
	    y_min = patch->getDomainLocalMin(1);
	    y_max = patch->getDomainLocalMax(1);
	}
	else {
	    y_min = max( y_min_global, patch->getDomainLocalMin(1) );
	    y_max = min( y_max_global, patch->getDomainLocalMax(1) );
	}
        if ( nDim_particle > 2 ) {
	    if (params.bc_em_type_trans=="periodic") {
		z_min = patch->getDomainLocalMin(2);
		z_max = patch->getDomainLocalMax(2);
	    }
	    else {
		z_min = max( z_min_global, patch->getDomainLocalMin(2) );
		z_max = min( z_max_global, patch->getDomainLocalMax(2) );
	    }
	}
    }

    // Define kind of boundary conditions
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
    else if ( params.species_param[ispec].bc_part_type_west == "adrien" ) {
	if (patch->isWestern()) bc_west = &adrien_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "none" ) {
	WARNING( "No Boundary Condition applied for species in west direction. Particles will be deleted. " << ispec );
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
    else if ( params.species_param[ispec].bc_part_type_east == "adrien" ) {
	if (patch->isEastern()) bc_east = &adrien_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "none" ) {
	WARNING( "No Boundary Condition applied for species in east direction. Particles will be deleted." << ispec );
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
	else if ( params.species_param[ispec].bc_part_type_south == "adrien" ) {
	    if (patch->isSouthern()) bc_south = &adrien_particle;
	}	
	else if ( params.species_param[ispec].bc_part_type_south == "none" ) {
	    WARNING( "No Boundary Condition applied for species in south direction. Particles will be deleted. " << ispec );
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
	else if ( params.species_param[ispec].bc_part_type_north == "adrien" ) {
	    if (patch->isNorthern()) bc_north = &adrien_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "none" ) {
	    WARNING( "No Boundary Condition applied for species in north direction. Particles will be deleted. " << ispec );
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
	}
    }

}

PartBoundCond::~PartBoundCond()
{
}

