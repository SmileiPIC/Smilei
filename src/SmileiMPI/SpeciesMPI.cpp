
#include <SpeciesMPI.h>

SpeciesMPI::SpeciesMPI()
{
}

void SpeciesMPI::init()
{
    int ndims_ = 2;
    int nbNeighbors_ = 2;

    for (int iDir=0 ; iDir<ndims_ ; iDir++) {
	for (int i=0 ; i<nbNeighbors_ ; i++) {
	    patch_buff_index_send[iDir][i].resize(0);
	    patch_buff_index_recv_sz[iDir][i] = 0;
	    patch_buff_index_send_sz[iDir][i] = 0;
	    //corner_buff_index_send[iDir][i].resize(0);
	    //corner_buff_index_recv_sz[iDir][i] = 0;
	    //corner_buff_index_send_sz[iDir][i] = 0;
	}
    }

    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    patchVectorRecv[iDim][iNeighbor].initialize(0,2);
	    patchVectorSend[iDim][iNeighbor].initialize(0,2);
	    // Globalize IP (store in SmileiMPI_Cart2D, OK while init+init+finalize / dir)
	    //cornerVectorRecv[iDim][iNeighbor].initialize(0,2);
	    //cornerVectorSend[iDim][iNeighbor].initialize(0,2);
	}
    }
}

SpeciesMPI::~SpeciesMPI()
{
}

