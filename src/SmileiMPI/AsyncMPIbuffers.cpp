
#include <AsyncMPIbuffers.h>

AsyncMPIbuffers::AsyncMPIbuffers()
{
}


AsyncMPIbuffers::~AsyncMPIbuffers()
{
}


void AsyncMPIbuffers::allocate(int ndims)
{
    srequest.resize(ndims);
    rrequest.resize(ndims);
    for (int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
    }
}

SpeciesMPIbuffers::SpeciesMPIbuffers()
{
}


SpeciesMPIbuffers::~SpeciesMPIbuffers()
{
}


void SpeciesMPIbuffers::allocate(int ndims)
{
    srequest.resize(ndims);
    rrequest.resize(ndims);

    partRecv.resize(ndims);
    partSend.resize(ndims);

    part_index_send.resize(ndims);
    part_index_send_sz.resize(ndims);
    part_index_recv_sz.resize(ndims);

    for (int i=0 ; i<ndims ; i++) {
        srequest[i].resize(2);
        rrequest[i].resize(2);
        partRecv[i].resize(2);
        partSend[i].resize(2);
        part_index_send[i].resize(2);
        part_index_send_sz[i].resize(2);
        part_index_recv_sz[i].resize(2);
    }

}

