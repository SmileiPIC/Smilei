
#include "AsyncMPIbuffers.h"
#include "Field.h"
#include "Patch.h"

#include <vector>
using namespace std;

AsyncMPIbuffers::AsyncMPIbuffers()
{
}


AsyncMPIbuffers::~AsyncMPIbuffers()
{
}


void AsyncMPIbuffers::allocate(unsigned int ndims)
{
    srequest.resize(ndims);
    rrequest.resize(ndims);
    for (unsigned int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
    }

    send_tags_.resize(ndims);
    recv_tags_.resize(ndims);
    for (unsigned int iDim = 0 ; iDim < ndims ; iDim++ ) {
        send_tags_[iDim].resize(2,MPI_PROC_NULL);
        recv_tags_[iDim].resize(2,MPI_PROC_NULL);
    }
}

void AsyncMPIbuffers::allocate(unsigned int ndims, Field* f, std::vector<unsigned int>& oversize)
{
    if (buf[0][0].size()!=0) return;
    srequest.resize(ndims);
    rrequest.resize(ndims);
    for (unsigned int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
    }
    
    std::vector<unsigned int> oversize2(oversize);
    oversize2[0] *= 2;
    oversize2[0] += 1 + f->isDual_[0];
    if (ndims>1) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f->isDual_[1];
        if (ndims>2) {
            oversize2[2] *= 2;
            oversize2[2] += 1 + f->isDual_[2];
        }
    }
    std::vector<unsigned int> n_elem = f->dims_;

    for (unsigned int iDim=0 ; iDim<ndims ; iDim++) {

        vector<int> idx( ndims,1 );
        idx[iDim] = 0;

        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {

            std::vector<unsigned int> tmp(3,1);
            tmp[0] =    idx[0]  * n_elem[0] + (1-idx[0]) * oversize2[0];
            if (ndims>1)
                tmp[1] =    idx[1]  * n_elem[1] + (1-idx[1]) * oversize2[1];
            if (ndims>2)
                tmp[2] =    idx[2]  * n_elem[2] + (1-idx[2]) * oversize2[2];

            buf[iDim][iNeighbor].resize( tmp[0]*tmp[1]*tmp[2], 0. );
        }
    }

    send_tags_.resize(ndims);
    recv_tags_.resize(ndims);
    for (unsigned int iDim = 0 ; iDim < ndims ; iDim++ ) {
        send_tags_[iDim].resize(2,MPI_PROC_NULL);
        recv_tags_[iDim].resize(2,MPI_PROC_NULL);
    }

}


void AsyncMPIbuffers::iallocate(unsigned int ndims, Field* f, std::vector<unsigned int>& oversize)
{
    if (ibuf[0][0].size()!=0) return;
    srequest.resize(ndims);
    rrequest.resize(ndims);
    for (unsigned int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
    }
    
    std::vector<unsigned int> oversize2(oversize);
    oversize2[0] *= 2;
    oversize2[0] += 1 + f->isDual_[0];
    if (ndims>1) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f->isDual_[1];
        if (ndims>2) {
            oversize2[2] *= 2;
            oversize2[2] += 1 + f->isDual_[2];
        }
    }
    std::vector<unsigned int> n_elem = f->dims_;

    for (unsigned int iDim=0 ; iDim<ndims ; iDim++) {

        vector<int> idx( ndims,1 );
        idx[iDim] = 0;

        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {

            std::vector<unsigned int> tmp(3,1);
            tmp[0] =    idx[0]  * n_elem[0] + (1-idx[0]) * oversize2[0];
            if (ndims>1)
                tmp[1] =    idx[1]  * n_elem[1] + (1-idx[1]) * oversize2[1];
            if (ndims>2)
                tmp[2] =    idx[2]  * n_elem[2] + (1-idx[2]) * oversize2[2];

            ibuf[iDim][iNeighbor].resize( tmp[0]*tmp[1]*tmp[2], 0. );
        }
    }

    send_tags_.resize(ndims);
    recv_tags_.resize(ndims);
    for (unsigned int iDim = 0 ; iDim < ndims ; iDim++ ) {
        send_tags_[iDim].resize(2,MPI_PROC_NULL);
        recv_tags_[iDim].resize(2,MPI_PROC_NULL);
    }

}


void AsyncMPIbuffers::defineTags(Patch* patch, int tag ) 
{
    for (unsigned int iDim=0 ; iDim< send_tags_.size() ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            send_tags_[iDim][iNeighbor] = buildtag( patch->hindex, iDim, iNeighbor, tag );
            recv_tags_[iDim][iNeighbor] = buildtag( patch->neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, tag );
        }

}


SpeciesMPIbuffers::SpeciesMPIbuffers()
{
}


SpeciesMPIbuffers::~SpeciesMPIbuffers()
{
}


void SpeciesMPIbuffers::allocate(unsigned int ndims)
{
    srequest.resize(ndims);
    rrequest.resize(ndims);

    partRecv.resize(ndims);
    partSend.resize(ndims);

    part_index_send.resize(ndims);
    part_index_send_sz.resize(ndims);
    part_index_recv_sz.resize(ndims);

    for (unsigned int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
        partRecv[i].resize(2);
        partSend[i].resize(2);
        part_index_send[i].resize(2);
        part_index_send_sz[i].resize(2);
        part_index_recv_sz[i].resize(2);
    }

}

