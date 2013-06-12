#include "Projector2D2Order.h"
#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particle.h"
#include "Tools.h"
#include "SmileiMPI_Cart2D.h"

#include <iostream>
#include <cmath>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Constructor for Projector2D2Order
// ---------------------------------------------------------------------------------------------------------------------
Projector2D2Order::Projector2D2Order (PicParams* params, SmileiMPI* smpi) : Projector2D(params, smpi)
{
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

	dx_inv_  = 1.0/params->cell_length[0];
	dx_ov_dt = params->cell_length[0] / params->timestep;

	index_domain_begin = smpi2D->getCellStartingGlobalIndex(0);

	DEBUG("cell_length "<< params->cell_length[0]);

}

Projector2D2Order::~Projector2D2Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd order projection in 1d3v simulations
// ---------------------------------------------------------------------------------------------------------------------
void Projector2D2Order::operator() (ElectroMagn* EMfields, Particle* part, double gf)
{
	MESSAGE( "to be implemented" );
    
} // END Projector2D2Order


void Projector2D2Order::operator() (Field* rho, Particle* part)
{
	MESSAGE( "to be implemented" );

}

