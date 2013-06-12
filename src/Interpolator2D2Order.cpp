
#include "Interpolator2D2Order.h"

#include <iostream>
#include <cmath>
using namespace std;

#include "ElectroMagn.h"
#include "Field2D.h"
#include "Particle.h" 
#include "SmileiMPI_Cart2D.h"

Interpolator2D2Order::Interpolator2D2Order(PicParams *params, SmileiMPI* smpi) : Interpolator2D(params, smpi)
{
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);

	dx_inv_ = 1.0/params->cell_length[0];

	index_domain_begin = smpi2D->getCellStartingGlobalIndex(0);

}

Interpolator2D2Order::~Interpolator2D2Order()
{
}

// ---------------------------------------------------------------------------------------------------------------------
// 2nd Order Interpolation of the fields at a the particle position (3 nodes are used)
// ---------------------------------------------------------------------------------------------------------------------
void Interpolator2D2Order::operator() (ElectroMagn* champs, Particle* part, LocalFields* ELoc, LocalFields* BLoc)
{
	MESSAGE( "to be implemented" );
	
} // END Interpolator2D2Order
