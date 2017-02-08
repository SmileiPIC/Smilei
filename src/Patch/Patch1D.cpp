
#include "Patch1D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch1D constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch1D::Patch1D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
  : Patch( params, smpi, ipatch, n_moved)
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCreation(params, smpi);
} // End Patch1D::Patch1D

// ---------------------------------------------------------------------------------------------------------------------
// Patch1D cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch1D::Patch1D(Patch1D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved, bool with_particles = true)
    : Patch( patch, params, smpi, ipatch, n_moved, with_particles )
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCloning(patch, params, smpi, with_particles);
} // End Patch1D::Patch1D


// ---------------------------------------------------------------------------------------------------------------------
// Patch1D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor 
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initStep2(Params& params)
{
    int xcall, ycall;
    
    Pcoordinates[0] = hindex;
    
    // 1st direction
    xcall = Pcoordinates[0]-1;
    ycall = Pcoordinates[1];
    if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
    neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    xcall = Pcoordinates[0]+1;
    if (params.bc_em_type_x[0]=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
    neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    
}

Patch1D::~Patch1D()
{
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        MPI_Type_free( &(ntype_[0][ix_isPrim]) );
        MPI_Type_free( &(ntype_[1][ix_isPrim]) );
        MPI_Type_free( &(ntypeSum_[0][ix_isPrim]) );
        }

}


void Patch1D::reallyinitSumField( Field* field, int iDim )
{
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initSumField( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(1, field);
    
    int tagp(0);
    if (field->name == "Jx") tagp = 1;
    if (field->name == "Jy") tagp = 2;
    if (field->name == "Jz") tagp = 3;
    if (field->name == "Rho") tagp = 4;

    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
    
    // Use a buffer per direction to exchange data before summing
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f1D->isDual_[0];
    
    int istart, ix;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
        
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]];
    
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
            ix = (1-iDim)*istart;
            int tag = buildtag( hindex, iDim, iNeighbor, tagp );
            MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.srequest[iDim][iNeighbor]) );
        } // END of Send
            
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            int tmp_elem = f1D->MPIbuff.buf[iDim][(iNeighbor+1)%2].size();
            int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, tagp );
            MPI_Irecv( &( f1D->MPIbuff.buf[iDim][(iNeighbor+1)%2][0]) , tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
        } // END of Recv
            
    } // END for iNeighbor
    
} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::finalizeSumField( Field* field, int iDim )
{
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f1D->isDual_[0];
    
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/

    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];
    
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &(f1D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            MPI_Wait( &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

        int istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
        int ix0 = (1-iDim)*istart;
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            int tmp = f1D->MPIbuff.buf[iDim][(iNeighbor+1)%2].size();
            for (unsigned int ix=0 ; ix< tmp ; ix++) {
                f1D->data_[ix0+ix] += f1D->MPIbuff.buf[iDim][(iNeighbor+1)%2][ix];
            }
        } // END if
    } // END for iNeighbor
}

void Patch1D::reallyfinalizeSumField( Field* field, int iDim )
{
} // END finalizeSumField


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initExchange( Field* field )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(1);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);

    int istart, ix;

    // Loop over dimField
    for (int iDim=0 ; iDim<nDim_fields_ ; iDim++) {

        MPI_Datatype ntype = ntype_[iDim][isDual[0]];
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
                ix = (1-iDim)*istart;
                int tag = send_tags_[iDim][iNeighbor];
                MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.srequest[iDim][iNeighbor]) );

            } // END of Send

            if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 - (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
                int tag = recv_tags_[iDim][iNeighbor];
                MPI_Irecv( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

            } // END of Recv

        } // END for iNeighbor

    } // END for iDim

} // END initExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI  (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::finalizeExchange( Field* field )
{
    Field1D* f1D =  static_cast<Field1D*>(field);

    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];

    // Loop over dimField
    for (int iDim=0 ; iDim<nDim_fields_ ; iDim++) {

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &(f1D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
             if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
               MPI_Wait( &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }

    } // END for iDim

} // END finalizeExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initExchange( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(1);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);

    int istart, ix;

    MPI_Datatype ntype = ntype_[iDim][isDual[0]];
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

            istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = (1-iDim)*istart;
            int tag = send_tags_[iDim][iNeighbor];
            MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.srequest[iDim][iNeighbor]) );

        } // END of Send

        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 - (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
            ix = (1-iDim)*istart;
            int tag = recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

        } // END of Recv

    } // END for iNeighbor

} // END initExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::finalizeExchange( Field* field, int iDim )
{
    Field1D* f1D =  static_cast<Field1D*>(field);

    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &(f1D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            MPI_Wait( &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }

} // END finalizeExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::createType( Params& params )
{
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual]
    int ny = oversize[0];
    int nline;

    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {

        // Standard Type
        ntype_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous(ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim]));    //line
        MPI_Type_commit( &(ntype_[0][ix_isPrim]) );

        ntype_[1][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous(clrw, MPI_DOUBLE, &(ntype_[1][ix_isPrim]));   //clrw lines
        MPI_Type_commit( &(ntype_[1][ix_isPrim]) );

        ntypeSum_[0][ix_isPrim] = MPI_DATATYPE_NULL;

        MPI_Datatype tmpType = MPI_DATATYPE_NULL;
        MPI_Type_contiguous(1, MPI_DOUBLE, &(tmpType));    //line
        MPI_Type_commit( &(tmpType) );


        nline = 1 + 2*params.oversize[0] + ix_isPrim;
        MPI_Type_contiguous(nline, tmpType, &(ntypeSum_[0][ix_isPrim]));    //line
        MPI_Type_commit( &(ntypeSum_[0][ix_isPrim]) );

        MPI_Type_free( &tmpType );
            
    }
    
} //END createType

