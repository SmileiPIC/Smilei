
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

}

