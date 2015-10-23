
#include <SpeciesMPI.h>

SpeciesMPI::SpeciesMPI()
{
}

void SpeciesMPI::init()
{
    // resize correctly arrays 
    // suppress corner object

    int ndims = 2;
    int nbNeighbors_ = 2;

    for (int iDir=0 ; iDir<ndims ; iDir++) {
	for (int i=0 ; i<nbNeighbors_ ; i++) {
	    patch_buff_index_send[iDir][i].resize(0);
	    patch_buff_index_recv_sz[iDir][i] = 0;
	    patch_buff_index_send_sz[iDir][i] = 0;
	    corner_buff_index_send[iDir][i].resize(0);
	    corner_buff_index_recv_sz[iDir][i] = 0;
	    corner_buff_index_send_sz[iDir][i] = 0;
	}
    }

    // ------------------------
    // Particles dimensions !!!
    // ------------------------
    for (int iDim=0 ; iDim<ndims ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    patchVectorRecv[iDim][iNeighbor].initialize(0,ndims);
	    patchVectorSend[iDim][iNeighbor].initialize(0,ndims);
	    // Globalize IP (store in SmileiMPI_Cart2D, OK while init+init+finalize / dir)
	    cornerVectorRecv[iDim][iNeighbor].initialize(0,ndims);
	    cornerVectorSend[iDim][iNeighbor].initialize(0,ndims);
	}
    }
}

SpeciesMPI::~SpeciesMPI()
{
}

