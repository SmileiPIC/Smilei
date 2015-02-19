#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( PicParams& params, int ispec, SmileiMPI* smpi )
{
    nDim_particle = params.nDim_particle;

    // Absolute global values
//    double x_min_global = params.cell_length[0]*n_ord_proj_max;
//    double x_max_global = params.cell_length[0]*( params.n_space_global[0]-1-n_ord_proj_max );
    double x_min_global = 0;
    double x_max_global = params.cell_length[0]*(params.n_space_global[0]);

    double y_min_global = 0;
    double y_max_global = params.cell_length[1]*(params.n_space_global[1]);
    double z_min_global = 0;
    double z_max_global = params.cell_length[2]*(params.n_space_global[2]);

    bc_west  = NULL;
    bc_east  = NULL;
    bc_south = NULL;
    bc_north = NULL;
    bc_bottom = NULL;
    bc_up     = NULL;

    // Define limits of local domain
    //if (!params.nspace_win_x) {
        x_min = max( x_min_global, smpi->getDomainLocalMin(0) );
        x_max = min( x_max_global, smpi->getDomainLocalMax(0) );
    //}
    //else {
    //    x_min = smpi->getDomainLocalMin(0);
    //    x_max = smpi->getDomainLocalMax(0);
    //}

    if ( nDim_particle > 1 ) {
	if (params.bc_em_type_trans=="periodic") {
	    y_min = smpi->getDomainLocalMin(1);
	    y_max = smpi->getDomainLocalMax(1);
	}
	else {
	    y_min = max( y_min_global, smpi->getDomainLocalMin(1) );
	    y_max = min( y_max_global, smpi->getDomainLocalMax(1) );
	}
        if ( nDim_particle > 2 ) {
	    if (params.bc_em_type_trans=="periodic") {
		z_min = smpi->getDomainLocalMin(2);
		z_max = smpi->getDomainLocalMax(2);
	    }
	    else {
		z_min = max( z_min_global, smpi->getDomainLocalMin(2) );
		z_max = min( z_max_global, smpi->getDomainLocalMax(2) );
	    }
	}
    }

    // Define kind of boundary conditions
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
    else if ( params.species_param[ispec].bc_part_type_west == "adrien" ) {
	if (smpi->isWestern()) bc_west = &adrien_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_west == "none" ) {
	WARNING( "No Boundary Condition applied for species in west direction " << ispec );
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
    else if ( params.species_param[ispec].bc_part_type_east == "adrien" ) {
	if (smpi->isEastern()) bc_east = &adrien_particle;
    }
    else if ( params.species_param[ispec].bc_part_type_east == "none" ) {
	WARNING( "No Boundary Condition applied for species in east direction " << ispec );
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
	else if ( params.species_param[ispec].bc_part_type_south == "adrien" ) {
	    if (smpi->isSouthern()) bc_south = &adrien_particle;
	}	
	else if ( params.species_param[ispec].bc_part_type_south == "none" ) {
	    WARNING( "No Boundary Condition applied for species in south direction " << ispec );
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
	else if ( params.species_param[ispec].bc_part_type_north == "adrien" ) {
	    if (smpi->isNorthern()) bc_north = &adrien_particle;
	}
	else if ( params.species_param[ispec].bc_part_type_north == "none" ) {
	    WARNING( "No Boundary Condition applied for species in north direction " << ispec );
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

void PartBoundCond::moveWindow_x(double shift, SmileiMPI* smpi)
{
    x_min += shift;
    x_max += shift;
    if (smpi->isWestern()) bc_west = &supp_particle;
    if (smpi->isEastern()) bc_east = &supp_particle;
}
