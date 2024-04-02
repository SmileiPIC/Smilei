
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


void AsyncMPIbuffers::allocate( unsigned int ndims )
{
    srequest.resize( ndims );
    rrequest.resize( ndims );
    for( unsigned int i=0 ; i<ndims ; i++ ) {
        srequest[i].resize( 2 );
        rrequest[i].resize( 2 );
    }
    
    send_tags_.resize( ndims );
    recv_tags_.resize( ndims );
    for( unsigned int iDim = 0 ; iDim < ndims ; iDim++ ) {
        send_tags_[iDim].resize( 2, MPI_PROC_NULL );
        recv_tags_[iDim].resize( 2, MPI_PROC_NULL );
    }
}

void AsyncMPIbuffers::defineTags( Patch *patch, SmileiMPI *smpi, int tag )
{
    for( unsigned int iDim=0 ; iDim< send_tags_.size() ; iDim++ )
        for( int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++ ) {
        
            int local_hindex = patch->hindex;
            if( patch->is_small ) {
                local_hindex -= smpi->patch_refHindexes[ patch->MPI_me_ ];
            }
            send_tags_[iDim][iNeighbor] = buildtag( local_hindex, iDim, iNeighbor, tag );
            
            local_hindex = patch->neighbor_[iDim][( iNeighbor+1 )%2];
            if( patch->is_small ) {
                if( patch->MPI_neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) {
                    local_hindex -= smpi->patch_refHindexes[ patch->MPI_neighbor_[iDim][( iNeighbor+1 )%2] ];
                }
            }
            if( patch->MPI_neighbor_[iDim][( iNeighbor+1 )%2]!=MPI_PROC_NULL ) {
                recv_tags_[iDim][iNeighbor] = buildtag( local_hindex, iDim, iNeighbor, tag );
            }
            
        }
        
}


SpeciesMPIbuffers::SpeciesMPIbuffers()
{
}


SpeciesMPIbuffers::~SpeciesMPIbuffers()
{
    for( size_t i=0 ; i<partRecv.size() ; i++ ) {
        delete partRecv[i][0];
        delete partRecv[i][1];
        delete partSend[i][0];
        delete partSend[i][1];
    }
}


void SpeciesMPIbuffers::allocate( unsigned int ndims )
{
    srequest.resize( ndims );
    rrequest.resize( ndims );
    
    partRecv.resize( ndims );
    partSend.resize( ndims );
    
    partSendSize.resize( ndims );
    partRecvSize.resize( ndims );
    
    for( unsigned int i=0 ; i<ndims ; i++ ) {
        srequest[i].resize( 2 );
        rrequest[i].resize( 2 );
        partRecv[i].resize( 2 );
        partSend[i].resize( 2 );
        partSendSize[i].resize( 2 );
        partRecvSize[i].resize( 2 );
        
        partRecv[i][0] = new Particles();
        partRecv[i][1] = new Particles();
        partSend[i][0] = new Particles();
        partSend[i][1] = new Particles();
    }
}

