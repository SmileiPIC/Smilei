
#include "AsyncMPIbuffers.h"
#include "ParticlesFactory.h"
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


void SpeciesMPIbuffers::allocate( Params &params, Patch *patch )
{
    srequest.resize( params.nDim_field );
    rrequest.resize( params.nDim_field );
    
    partRecv.resize( params.nDim_field );
    partSend.resize( params.nDim_field );
    
    partSendSize.resize( params.nDim_field );
    partRecvSize.resize( params.nDim_field );
    
    for( unsigned int i=0 ; i<params.nDim_field ; i++ ) {
        srequest[i].resize( 2 );
        rrequest[i].resize( 2 );
        partRecvSize[i].resize( 2 );
        partSendSize[i].resize( 2 );
        
        // NOTE: send/recv buffers on xmin / xmax use a different constructor because
        //       they must be sent on GPU for exchanging particles
        partRecv[i].resize( 2 );
        partSend[i].resize( 2 );
        if( i == 0 ) {
            partRecv[i][0] = ParticlesFactory::create( params, *patch );
            partRecv[i][1] = ParticlesFactory::create( params, *patch );
            partSend[i][0] = ParticlesFactory::create( params, *patch );
            partSend[i][1] = ParticlesFactory::create( params, *patch );
        } else {
            partRecv[i][0] = new Particles();
            partRecv[i][1] = new Particles();
            partSend[i][0] = new Particles();
            partSend[i][1] = new Particles();
        }
    }
}

