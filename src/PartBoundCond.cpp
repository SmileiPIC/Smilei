#include "PartBoundCond.h"

#include <cstdlib>

#include <iostream>
#include <string>

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "SmileiMPI.h"
#include "Tools.h"

using namespace std;

PartBoundCond::PartBoundCond( PicParams *params, int ispec, SmileiMPI* smpi )
{
    nDim_particle = params->nDim_particle;

    int n_ord_proj_max = 5;
    //!\todo define n_ord_proj_max, value from SQUASH
    if (params->interpolation_order==2) n_ord_proj_max = 5;
    if (params->interpolation_order==3) n_ord_proj_max = 5;

    // Absolute global values
//    double x_min_global = params->cell_length[0]*n_ord_proj_max;
//    double x_max_global = params->cell_length[0]*( params->n_space_global[0]-1-n_ord_proj_max );
    double x_min_global = 0;
    double x_max_global = params->cell_length[0]*(params->n_space_global[0]);

    double y_min_global = 0;
    double y_max_global = params->cell_length[1]*(params->n_space_global[1]);
    double z_min_global = 0;
    double z_max_global = params->cell_length[2]*(params->n_space_global[2]);

    bc_west  = NULL;
    bc_east  = NULL;
    bc_south = NULL;
    bc_north = NULL;
    bc_bottom = NULL;
    bc_up     = NULL;

    // Define limits of local domain
    //if (!params->res_space_win_x) {
        x_min = max( x_min_global, smpi->getDomainLocalMin(0) );
        x_max = min( x_max_global, smpi->getDomainLocalMax(0) );
    //}
    //else {
    //    x_min = smpi->getDomainLocalMin(0);
    //    x_max = smpi->getDomainLocalMax(0);
    //}

    if ( nDim_particle > 1 ) {
	if (params->bc_em_type_trans=="periodic") {
	    y_min = smpi->getDomainLocalMin(1);
	    y_max = smpi->getDomainLocalMax(1);
	}
	else {
	    y_min = max( y_min_global, smpi->getDomainLocalMin(1) );
	    y_max = min( y_max_global, smpi->getDomainLocalMax(1) );
	}
        if ( nDim_particle > 2 ) {
	    if (params->bc_em_type_trans=="periodic") {
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
    if ( params->species_param[ispec].bc_part_type_long == "refl" ) {
        if (smpi->isWestern()) bc_west = &refl_particle;
        if (smpi->isEastern()) bc_east = &refl_particle;
    }
    else if ( params->species_param[ispec].bc_part_type_long == "supp" ) {
        if (smpi->isWestern()) bc_west = &supp_particle;
        if (smpi->isEastern()) bc_east = &supp_particle;
    }
    else if ( params->species_param[ispec].bc_part_type_long == "stop" ) {
        if (smpi->isWestern()) bc_west = &stop_particle;
        if (smpi->isEastern()) bc_east = &stop_particle;
    }
    else if ( params->species_param[ispec].bc_part_type_long == "none" ) {
        WARNING( "No Boundary Condition applied for species in longitudinal direction " << ispec );
    }
    else {
	ERROR( "Longitudinal boundary condition undefined" );
    }

    if ( nDim_particle > 1 ) {
	//if  (params->bc_em_type_trans!="periodic") {
	if ( params->species_param[ispec].bc_part_type_trans == "refl" ) {
	    if (smpi->isSouthern()) bc_south = &refl_particle;
	    if (smpi->isNorthern()) bc_north = &refl_particle;
	}
	else if ( params->species_param[ispec].bc_part_type_trans == "supp" ) {
	    if (smpi->isSouthern()) bc_south = &supp_particle;
	    if (smpi->isNorthern()) bc_north = &supp_particle;
	}
	else if ( params->species_param[ispec].bc_part_type_trans == "stop" ) {
	    if (smpi->isSouthern()) bc_south = &stop_particle;
	    if (smpi->isNorthern()) bc_north = &stop_particle;
	}
	else if ( params->species_param[ispec].bc_part_type_trans == "none" ) {
	    WARNING( "No Boundary Condition applied for species in transverse direction " << ispec );
	}
	else {
	    ERROR( "Transverse boundary condition undefined : " << params->species_param[ispec].bc_part_type_trans  );
	}
	//} // else NULL
	if ( nDim_particle > 2 ) {
	    //if  (params->bc_em_type_trans!="periodic") {
	    if ( params->species_param[ispec].bc_part_type_trans == "refl" ) {
		if (z_min==z_min_global) bc_bottom = &refl_particle;
		if (z_max==z_max_global) bc_up     = &refl_particle;
	    }
	    else if ( params->species_param[ispec].bc_part_type_trans == "supp" ) {
		if (z_min==z_min_global) bc_bottom = &supp_particle;
		if (z_max==z_max_global) bc_up     = &supp_particle;
	    }
	    else if ( params->species_param[ispec].bc_part_type_trans == "stop" ) {
		if (z_min==z_min_global) bc_bottom = &stop_particle;
		if (z_max==z_max_global) bc_up     = &stop_particle;
	    }
	    //} // else NULL
	}
    }

}

PartBoundCond::~PartBoundCond()
{
}

void PartBoundCond::moveWindow_x(double shift)
{
    x_min += shift;
    x_max += shift;

}
