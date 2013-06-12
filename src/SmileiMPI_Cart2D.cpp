#include "SmileiMPI_Cart2D.h"

#include "Species.h"
#include "ParticleFactory.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include "Tools.h" 

#include <string>
using namespace std;

SmileiMPI_Cart2D::SmileiMPI_Cart2D( int* argc, char*** argv )
	: SmileiMPI( argc, argv )
{
}

SmileiMPI_Cart2D::SmileiMPI_Cart2D( SmileiMPI* smpi)
	: SmileiMPI( smpi )
{
	ndims_ = 2;
	dims_  = new int(ndims_);
	coords_  = new int(ndims_);
	periods_  = new int(ndims_);
	reorder_ = 0;

	nbNeighbors_ = 2;
	neighbor_  = new int(nbNeighbors_);

	buff_send = new std::vector<Particle*>[nbNeighbors_];
	buff_recv = new std::vector<Particle*>[nbNeighbors_];

	for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
	for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
	for (int i=0 ; i<ndims_ ; i++) dims_[i] = 0;

	for (int i=0 ; i<nbNeighbors_ ; i++) {
		neighbor_[i] = MPI_PROC_NULL;
		buff_send[i].resize(0);
		buff_recv[i].resize(0);
	}
}

SmileiMPI_Cart2D::~SmileiMPI_Cart2D()
{
	delete dims_;
	delete periods_;
	delete coords_;
	delete neighbor_;

	delete [] buff_send;
	delete [] buff_recv;

	if ( SMILEI_COMM_2D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_2D);

}

void SmileiMPI_Cart2D::createTopology()
{
	MPI_Dims_create( smilei_sz, ndims_, dims_ );
	MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, dims_, periods_, reorder_, &SMILEI_COMM_2D );
	MPI_Cart_coords( SMILEI_COMM_2D, smilei_rk, ndims_, coords_ );

	// neighbor_[0]  |  Current process  |  neighbor_[1] //
	MPI_Cart_shift( SMILEI_COMM_2D, 0, 1, &(neighbor_[0]), &(neighbor_[1]) );
	PMESSAGE ( 0, smilei_rk, "Neighbors of process : " << neighbor_[0] << " - " << neighbor_[1]  );

}

void SmileiMPI_Cart2D::exchangeParticles(Species* species, int ispec, PicParams* params)
{
	MESSAGE( "to be implemented" );

} // END exchangeParticles


void SmileiMPI_Cart2D::sumField( Field* field )
{
	MESSAGE( "to be implemented" );


} // END sumField


void SmileiMPI_Cart2D::exchangeField( Field* field )
{
	MESSAGE( "to be implemented" );


} // END exchangeField


void SmileiMPI_Cart2D::writeField( Field* field, string name )
{
	MESSAGE( "to be implemented" );

} // END writeField
