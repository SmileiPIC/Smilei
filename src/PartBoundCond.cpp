
#include "PartBoundCond.h"

#include "Particle.h" 
#include "BoundaryConditionType.h" 
#include "Tools.h" 

#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

PartBoundCond::PartBoundCond( PicParams *params, int ispec )
{
	//!\todo deine n_ord_proj_max 
        int n_ord_proj_max = 5;

	x_min = 0.;
	x_max = 0.;
	y_min = 0.;
	y_max = 0.;
	z_min = 0.;
	z_max = 0.;

	bc_east  = NULL;
	bc_west  = NULL;
	bc_north  = NULL;
	bc_south  = NULL;
	bc_bottom = NULL;
	bc_up     = NULL;
	
	if ( params->nDim_particle == 1 ) {
        	x_min = params->cell_length[0]*n_ord_proj_max;
	        x_max = params->cell_length[0]*( params->n_space[ispec]+1-n_ord_proj_max );
		if ( params->species_param[ispec].bc_part_type == "refl" ) {
			bc_east = &refl_particle;
			bc_west = &refl_particle;
		}
		else {
			WARNING( "No Boundary Condition applied for species " << ispec );
		}
	}
	else if ( params->nDim_particle == 2 ) {
		WARNING ( "Allocation of 2D BC particules not yet implemented !" ); 
	}
	else if ( params->nDim_particle == 3 ) {
		WARNING ( "Allocation of 3D BC particules not yet implemented !" ); 
	}
	else {
		ERROR( "Bad Particle Dimensions" );
	} 
	
}

PartBoundCond::~PartBoundCond()
{
}

