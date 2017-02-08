
#include <AsyncMPIbuffers.h>
#include <Field.h>

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
}

void AsyncMPIbuffers::allocate(int ndims, Field* f)
{
    if (srequest.size()!=0) return;
    srequest.resize(ndims);
    rrequest.resize(ndims);
    for (unsigned int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
    }
    
    std::vector<unsigned int> oversize2(ndims,2);
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

    for (int iDim=0 ; iDim<ndims ; iDim++) {

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

