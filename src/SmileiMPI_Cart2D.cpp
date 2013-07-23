#include "SmileiMPI_Cart2D.h"

#include "Species.h"
#include "ParticleFactory.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include "Tools.h" 

#include <string>
#include <mpi.h>
#include <cmath>

using namespace std;

SmileiMPI_Cart2D::SmileiMPI_Cart2D( int* argc, char*** argv )
	: SmileiMPI( argc, argv )
{
}

SmileiMPI_Cart2D::SmileiMPI_Cart2D( SmileiMPI* smpi)
	: SmileiMPI( smpi )
{
	ndims_ = 2;
	number_of_procs = new int(ndims_);
	coords_  = new int(ndims_);
	periods_  = new int(ndims_);
	reorder_ = 0;

	nbNeighbors_ = 2; // per direction

	for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
	// Geometry periodic in y
	periods_[1] = 1;
	if (periods_[1] == 1)
		PMESSAGE( 0, smilei_rk, "Periodic geometry / y");
	for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
	for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;

	for (int iDim=0 ; iDim<ndims_ ; iDim++)
		for (int iNeighbors=0 ; iNeighbors<nbNeighbors_ ; iNeighbors++)
			neighbor_[iDim][iNeighbors] = MPI_PROC_NULL;

	for (int iDim=0 ; iDim<ndims_ ; iDim++) {
		for (int i=0 ; i<nbNeighbors_ ; i++) {
			buff_index_send[iDim][i].resize(0);
			buff_index_recv_sz[iDim][i] = 0;
		}
	}

}

SmileiMPI_Cart2D::~SmileiMPI_Cart2D()
{
	for (int ix_isPrim=0 ; ix_isPrim<1 ; ix_isPrim++) {
		for (int iy_isPrim=0 ; iy_isPrim<1 ; iy_isPrim++) {
			MPI_Type_free( &ntype_[0][ix_isPrim][iy_isPrim]); //line
			MPI_Type_free( &ntype_[1][ix_isPrim][iy_isPrim]); // column

			MPI_Type_free( &ntypeSum_[0][ix_isPrim][iy_isPrim]); //line
			MPI_Type_free( &ntypeSum_[1][ix_isPrim][iy_isPrim]); // column
		}
	}

	delete number_of_procs;
	delete periods_;
	delete coords_;

	if ( SMILEI_COMM_2D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_2D);

}

void SmileiMPI_Cart2D::createTopology(PicParams& params)
{
	if (params.nDim_field == 2) {
		double tmp = params.res_space[0]*params.sim_length[0] / ( params.res_space[1]*params.sim_length[1] );
		number_of_procs[0] = min( smilei_sz, max(1, (int)sqrt ( (double)smilei_sz*tmp*tmp) ) );
		number_of_procs[1] = (int)(smilei_sz / number_of_procs[0]);
	}
//	number_of_procs[0] = 3;
//	number_of_procs[1] = 1;

	MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_2D );
	MPI_Cart_coords( SMILEI_COMM_2D, smilei_rk, ndims_, coords_ );

	//                  |                   |                  //
	//                  |  neighbor_[2][1]  |                  //
	//                  |                   |                  //

	//                  |  neighbor_[1][1]  |                  //
	// neighbor_[0][0]  |  Current process  |  neighbor_[0][1] //
	//                  |  neighbor_[1][0]  |                  //

	//                  |                   |                  //
	//                  |  neighbor_[2][0]  |                  //
	//                  |                   |                  //

	// ==========================================================
	// ==========================================================
	// ==========================================================

	// crossNei_[x][x]  | crossNei_[x][x]   | crossNei_[x][x]  //
	// crossNei_[x][x]  |                   | crossNei_[x][x]  //
	// crossNei_[x][x]  | crossNei_[x][x]   | crossNei_[x][x]  //

	// crossNei_[x][x]  |                   | crossNei_[x][x]  //
	//                  |  Current process  |                  //		-> Manage working direction per direction
	// crossNei_[x][x]  |                   | crossNei_[x][x]  //

	// crossNei_[x][x]  | crossNei_[x][x]   | crossNei_[x][x]  //
	// crossNei_[x][x]  |                   | crossNei_[x][x]  //
	// crossNei_[x][x]  | crossNei_[x][x]   | crossNei_[x][x]  //


	for (int iDim=0 ; iDim<ndims_ ; iDim++) {
		MPI_Cart_shift( SMILEI_COMM_2D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
		PMESSAGE ( 0, smilei_rk, "Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " - " << neighbor_[iDim][1]  );
	}


	for (unsigned int i=0 ; i<params.nDim_field ; i++) {

		params.n_space[i] = params.n_space_global[i] / number_of_procs[i];

		n_space_global[i] = params.n_space_global[i];
		oversize[i] = params.oversize[i] = 2;
		cell_starting_global_index[i] = coords_[i]*(params.n_space_global[i] / number_of_procs[i]);
		// min/max_local : describe local domain in which particles cat be moved
		//                 different from domain on which E, B, J are defined
		min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
		max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
		cell_starting_global_index[i] -= params.oversize[i];

		if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
			// Correction on the last MPI process of the direction to use the wished number of cells
			if (coords_[i]==number_of_procs[i]-1) {
				params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
			}
		}

	}
	MESSAGE( "n_space / rank " << smilei_rk << " = " << params.n_space[0] << " " << params.n_space[1] );


}

void SmileiMPI_Cart2D::exchangeParticles(Species* species, int ispec, PicParams* params)
{
	std::vector<Particle*>* cuParticles = &species->particles;

	/********************************************************************************/
	// Build lists of indexes of particle to exchange per neighbor
	// Computed from indexes_of_particles_to_exchange computed during particles' BC
	/********************************************************************************/
	int n_part_send = indexes_of_particles_to_exchange.size();
	int n_part_recv;

	int iPart;
	for (int i=0 ; i<n_part_send ; i++) {
		iPart = indexes_of_particles_to_exchange[i];
		for (int iDim=0 ; iDim<ndims_ ; iDim++) {
			if ( (*cuParticles)[iPart]->position(iDim) < min_local[iDim]) {
				buff_index_send[iDim][0].push_back( indexes_of_particles_to_exchange[i] );
				break;
			}
			if ( (*cuParticles)[iPart]->position(iDim) >= max_local[iDim]) {
				buff_index_send[iDim][1].push_back( indexes_of_particles_to_exchange[i] );
				break;
			}
		}
	} // END for iPart = f(i)


	/********************************************************************************/
	// Exchange particles
	/********************************************************************************/
	for (int iDim=0 ; iDim<ndims_ ; iDim++) {

		MPI_Status stat[2];
		MPI_Request request[2];

		/********************************************************************************/
		// Exchange number of particles to exchange to establish or not a communication
		/********************************************************************************/
		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
				n_part_send = (buff_index_send[iDim][iNeighbor]).size();
				MPI_Isend( &n_part_send, 1, MPI_INT, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(request[iNeighbor]) );
			} // END of Send
			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				MPI_Irecv( &(buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(request[(iNeighbor+1)%2]) );
			}
		}

		/********************************************************************************/
		// Wait for end of communications over number of particles
		/********************************************************************************/
		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[iNeighbor]), &(stat[iNeighbor]) );
			}
			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[(iNeighbor+1)%2]), &(stat[(iNeighbor+1)%2]) );
			}
		}


		/********************************************************************************/
		// Define buffers to exchange buff_index_send[iDim][iNeighbor].size();
		/********************************************************************************/
		std::vector<unsigned int> partSize(2,0);
		partSize[0] = 0; // Number of particle exchanged per direction
		//! \todo{6 replaced by number of properties by particle}
		partSize[1] = 6;
		Field2D partArrayRecv[2];
		Field2D partArraySend[2];

		/********************************************************************************/
		// Proceed to effective Particles' communications
		/********************************************************************************/
		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

			// n_part_send : number of particles to send to current neighbor
			n_part_send = (buff_index_send[iDim][iNeighbor]).size();
			if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
				partSize[0] = n_part_send;
				partArraySend[iNeighbor].allocateDims(partSize);
				for (int iPart=0 ; iPart<n_part_send ; iPart++) {
					memcpy(&(partArraySend[iNeighbor](iPart,0)), &((*cuParticles)[ buff_index_send[iDim][iNeighbor][iPart] ]->position(0)), 6*sizeof(double) );
				}
				MPI_Isend( &(partArraySend[iNeighbor](0,0)), (int)6*n_part_send, MPI_DOUBLE, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(request[iNeighbor]) );

			} // END of Send

			n_part_recv = buff_index_recv_sz[iDim][(iNeighbor+1)%2];
			if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
				partSize[0] = n_part_recv;
				partArrayRecv[(iNeighbor+1)%2].allocateDims(partSize);
				MPI_Irecv( &(partArrayRecv[(iNeighbor+1)%2](0,0)), (int)6*n_part_recv, MPI_DOUBLE,  neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(request[(iNeighbor+1)%2]) );
			} // END of Recv

		} // END for iNeighbor


		/********************************************************************************/
		// Wait for end of communications over Particles
		/********************************************************************************/
		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

			n_part_send = buff_index_send[iDim][iNeighbor].size();
			if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ){
				MPI_Wait( &(request[iNeighbor]), &(stat[iNeighbor]) );
			}

			n_part_recv = buff_index_recv_sz[iDim][(iNeighbor+1)%2];
			if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
				MPI_Wait( &(request[(iNeighbor+1)%2]), &(stat[(iNeighbor+1)%2]) );
				int n_particles = species->getNbrOfParticles();
				cuParticles->resize( n_particles + n_part_recv );
				for (int iPart=0 ; iPart<n_part_recv; iPart++ ) {
					(*cuParticles)[n_particles+iPart] = ParticleFactory::create(params, ispec);
					memcpy( &( ((*cuParticles)[n_particles+iPart])->position(0) ), &(partArrayRecv[(iNeighbor+1)%2](iPart,0)),6*sizeof(double) );
				}
			}

		}

		/********************************************************************************/
		// Clean lists of indexes of particle to exchange per neighbor
		/********************************************************************************/
		for (int i=0 ; i<nbNeighbors_ ; i++)
			buff_index_send[iDim][i].clear();

	}

	/********************************************************************************/
	// Delete Particles included in buff_send/buff_recv
	/********************************************************************************/
	for (int i=n_part_send-1 ; i>=0 ; i--) {
		iPart = indexes_of_particles_to_exchange[i];
		(*cuParticles)[iPart]->~Particle();
		cuParticles->erase( cuParticles->begin() + iPart );
	} // END for iPart = f(i)

} // END exchangeParticles


void SmileiMPI_Cart2D::createType( PicParams& params )
{
	int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
	int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];

	// MPI_Datatype ntype_[nDim][primDual][primDual]
	int nx, ny;
	int nline, ncol;
	for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
		nx = nx0 + ix_isPrim;
		for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
			ny = ny0 + iy_isPrim;
			ntype_[0][ix_isPrim][iy_isPrim] = NULL;
			MPI_Type_contiguous(ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim][iy_isPrim]));    //line
			MPI_Type_commit( &(ntype_[0][ix_isPrim][iy_isPrim]) );
			ntype_[1][ix_isPrim][iy_isPrim] = NULL;
			MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &(ntype_[1][ix_isPrim][iy_isPrim])); // column
			MPI_Type_commit( &(ntype_[1][ix_isPrim][iy_isPrim]) );

			ntypeSum_[0][ix_isPrim][iy_isPrim] = NULL;
			nline = 1 + 2*params.oversize[0] + ix_isPrim;
			MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
			MPI_Type_commit( &(ntypeSum_[0][ix_isPrim][iy_isPrim]) );
			ntypeSum_[1][ix_isPrim][iy_isPrim] = NULL;
			ncol  = 1 + 2*params.oversize[1] + iy_isPrim;
			MPI_Type_vector(nx, ncol, ny, MPI_DOUBLE, &(ntypeSum_[1][ix_isPrim][iy_isPrim])); // column
			MPI_Type_commit( &(ntypeSum_[1][ix_isPrim][iy_isPrim]) );

		}
	}

}


void SmileiMPI_Cart2D::sumField( Field* field )
{
	std::vector<unsigned int> n_elem = field->dims_;
	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field2D* f2D =  static_cast<Field2D*>(field);


	// Use a buffer per direction to exchange data before summing
	Field2D buf[ndims_][ nbNeighbors_ ];
	// Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
	std::vector<unsigned int> oversize2 = oversize;
	oversize2[0] *= 2; oversize2[0] += 1 + f2D->isPrimal_[0];
	oversize2[1] *= 2; oversize2[1] += 1 + f2D->isPrimal_[1];

	for (int iDim=0 ; iDim<ndims_ ; iDim++) {
		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			std::vector<unsigned int> tmp(ndims_,0);
			tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
			tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
			buf[iDim][iNeighbor].allocateDims( tmp );
		}
	}

	int istart, ix, iy;

	/********************************************************************************/
	// Send/Recv in a buffer data to sum
	/********************************************************************************/
	for (int iDim=0 ; iDim<ndims_ ; iDim++) {

		MPI_Datatype ntype = ntypeSum_[iDim][isPrimal[0]][isPrimal[1]];
		MPI_Status stat[2];
		MPI_Request request[2];

		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
				istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
				ix = (1-iDim)*istart;
				iy =    iDim *istart;
				MPI_Isend( &(f2D->data_[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(request[iNeighbor]) );
			} // END of Send

			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
				MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_[0][0] ), tmp_elem, MPI_DOUBLE, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(request[(iNeighbor+1)%2]) );
			} // END of Recv

		} // END for iNeighbor


		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[iNeighbor]), &(stat[iNeighbor]) );
			}
			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[(iNeighbor+1)%2]), &(stat[(iNeighbor+1)%2]) );
			}
		}


		// Synchro before summing, to not sum with data ever sum
		// Merge loops, Sum direction by direction permits to not communicate with diagonal neighbors
		barrier();
		/********************************************************************************/
		// Sum data on each process, same operation on both side
		/********************************************************************************/

		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
			int ix0 = (1-iDim)*istart;
			int iy0 =    iDim *istart;
			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
					for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
						f2D->data_[ix0+ix][iy0+iy] += (buf[iDim][(iNeighbor+1)%2])(ix,iy);
				}
			} // END if

		} // END for iNeighbor
		
		barrier();
		
	} // END for iDim

} // END sumField


void SmileiMPI_Cart2D::exchangeField( Field* field )
{
	std::vector<unsigned int> n_elem   = field->dims_;
	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field2D* f2D =  static_cast<Field2D*>(field);

	int istart, ix, iy;

	// Loop over dimField
	for (int iDim=0 ; iDim<ndims_ ; iDim++) {

		MPI_Datatype ntype = ntype_[iDim][isPrimal[0]][isPrimal[1]];
		MPI_Status stat[2];
		MPI_Request request[2];

		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {

				istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isPrimal[iDim]) ) + (1-iNeighbor) * ( 2*oversize[iDim]+1-(1-isPrimal[iDim]) );
				ix = (1-iDim)*istart;
				iy =    iDim *istart;
				MPI_Isend( &(f2D->data_[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(request[iNeighbor]) );

			} // END of Send

			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {

				istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
				ix = (1-iDim)*istart;
				iy =    iDim *istart;
				MPI_Irecv( &(f2D->data_[ix][iy]), 1, ntype, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(request[(iNeighbor+1)%2]));

			} // END of Recv

		} // END for iNeighbor

		for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
			if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[iNeighbor]), &(stat[iNeighbor]) );
			}
			if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
				MPI_Wait( &(request[(iNeighbor+1)%2]), &(stat[(iNeighbor+1)%2]) );
			}
		}

	} // END for iDim


} // END exchangeField
