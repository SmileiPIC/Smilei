
#include "PartBoundCond.h"

#include "Particles.h"
#include "BoundaryConditionType.h"
#include "SmileiMPI.h"
#include "Tools.h"

#include <iostream>
#include <string>
#include <cstdlib>
using namespace std;

PartBoundCond::PartBoundCond( PicParams *params, int ispec, SmileiMPI* smpi )
{
    nDim_particle = params->nDim_particle;

    int n_ord_proj_max = 5;
    //!\todo define n_ord_proj_max, value from SQUASH
    if (params->interpolation_order==2) n_ord_proj_max = 5;
    if (params->interpolation_order==3) n_ord_proj_max = 5;

    // Absolute global values
    double x_min_global = params->cell_length[0]*n_ord_proj_max;
    double x_max_global = params->cell_length[0]*( params->n_space_global[0]-1-n_ord_proj_max );
    double y_min_global = params->cell_length[1]*n_ord_proj_max;
    double y_max_global = params->cell_length[1]*( params->n_space_global[1]-1-n_ord_proj_max );
    double z_min_global = params->cell_length[2]*n_ord_proj_max;
    double z_max_global = params->cell_length[2]*( params->n_space_global[2]-1-n_ord_proj_max );

    bc_west  = NULL;
    bc_east  = NULL;
    bc_south = NULL;
    bc_north = NULL;
    bc_bottom = NULL;
    bc_up     = NULL;

    x_min = max( x_min_global, smpi->getDomainLocalMin(0) );
    x_max = min( x_max_global, smpi->getDomainLocalMax(0) );
    if ( nDim_particle > 1 ) {
        // if ( smpi->periods_[1] == 1) {
        y_min = smpi->getDomainLocalMin(1);
        y_max = smpi->getDomainLocalMax(1);
//		}
//		else {
//			y_min = max( y_min_global, smpi->getDomainLocalMin(1) );
//			y_max = min( y_max_global, smpi->getDomainLocalMax(1) );
//		}
        if ( nDim_particle > 2 ) {
            z_min = max( z_min_global, smpi->getDomainLocalMin(2) );
            z_max = min( z_max_global, smpi->getDomainLocalMax(2) );
        }
    }

    if ( params->species_param[ispec].bc_part_type == "refl" ) {
        if (x_min==x_min_global) bc_west = &refl_particle;
        if (x_max==x_max_global) bc_east = &refl_particle;
        if ( nDim_particle > 1 ) {
//			if ( smpi->periods_[1] != 1) {
//				if (y_min==y_min_global) bc_south = &refl_particle;
//				if (y_max==y_max_global) bc_north = &refl_particle;
//			}
            if ( nDim_particle > 2 ) {
                if (z_min==z_min_global) bc_bottom = &refl_particle;
                if (z_max==z_max_global) bc_up     = &refl_particle;
            }
        }
    }
    else
        WARNING( "No Boundary Condition applied for species " << ispec );

}

PartBoundCond::~PartBoundCond()
{
}

