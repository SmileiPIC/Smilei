
#include "Patch1D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"
#include "SmileiIOFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch1D constructor :
//   - Pcoordinates, neighbor_ resized in Patch constructor 
//   - Call Patch::finalizePatchInit to allocate data structure
// ---------------------------------------------------------------------------------------------------------------------
Patch1D::Patch1D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
  : Patch( params, smpi, ipatch, n_moved)
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

    // Call generic Patch::finalizePatchInit method
    finalizePatchInit( params, smpi, n_moved );

} // End Patch1D::Patch1D


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initSumField( Field* field, int iDim )
{
    vector<unsigned int> patch_oversize(nDim_fields_,2);
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = patch_oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f1D->isDual_[0];
    
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	std::vector<unsigned int> tmp(nDim_fields_,0);
	tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
	buf[iDim][iNeighbor].allocateDims( tmp );
    }
     
    int istart, ix;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
        
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]];
        
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            
	if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
	    istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
	    ix = (1-iDim)*istart;
	    int tag = buildtag( hindex, iDim, iNeighbor );
	    MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_srequest[iDim][iNeighbor]) );
	    //MPI_Isend( &(f1D->data_[ix]), iDim  * n_elem[0] + (1-iDim) * oversize2[0], MPI_DOUBLE, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_srequest[iDim][iNeighbor]) );
	} // END of Send
            
	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
	    int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0];
	    int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
	    MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_[0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
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
    vector<unsigned int> patch_oversize(nDim_fields_,2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = patch_oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f1D->isDual_[0];
    
    int istart, ix;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/

    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]];
    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];
	
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
	    MPI_Wait( &(f1D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
	}
	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
	    MPI_Wait( &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
	}
    }


    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
       
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
	int ix0 = (1-iDim)*istart;
	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
	    for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
		f1D->data_[ix0+ix] += (buf[iDim][(iNeighbor+1)%2])(ix);
	    }
	} // END if
            
    } // END for iNeighbor
        

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	buf[iDim][iNeighbor].deallocateDims();
    }

} // END finalizeSumField


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initExchange( Field* field )
{
    vector<unsigned int> patch_oversize(nDim_fields_,2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);

    int istart, ix;

    // Loop over dimField
    for (int iDim=0 ; iDim<nDim_fields_ ; iDim++) {

        MPI_Datatype ntype = ntype_[iDim][isDual[0]];
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

	    if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*patch_oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*patch_oversize[iDim] + isDual[iDim] );
                ix = (1-iDim)*istart;
		int tag = buildtag( hindex, iDim, iNeighbor );
                MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_srequest[iDim][iNeighbor]) );

            } // END of Send

	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
 		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
		MPI_Irecv( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));

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
                MPI_Wait( &(f1D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
 	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
               MPI_Wait( &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
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
    vector<unsigned int> patch_oversize(nDim_fields_,2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);

    int istart, ix;

    MPI_Datatype ntype = ntype_[iDim][isDual[0]];
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

	if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

	    istart = iNeighbor * ( n_elem[iDim]- (2*patch_oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*patch_oversize[iDim] + isDual[iDim] );
	    ix = (1-iDim)*istart;
	    int tag = buildtag( hindex, iDim, iNeighbor );
	    MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_srequest[iDim][iNeighbor]) );

	} // END of Send

	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

	    istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
	    ix = (1-iDim)*istart;
	    int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
	    MPI_Irecv( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));

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
	    MPI_Wait( &(f1D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
	}
	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
	    MPI_Wait( &(f1D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
	}
    }

} // END finalizeExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::createType( Params& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual]
    int nx;
    int ny = 1;
    int nline;

    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;

	// Standard Type
	ntype_[0][ix_isPrim] = NULL;
	MPI_Type_contiguous(ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim]));    //line
	MPI_Type_commit( &(ntype_[0][ix_isPrim]) );

	ntype_[1][ix_isPrim] = NULL;
	MPI_Type_contiguous(ny*clrw, MPI_DOUBLE, &(ntype_[1][ix_isPrim]));   //clrw lines
	MPI_Type_commit( &(ntype_[1][ix_isPrim]) );

	ntypeSum_[0][ix_isPrim] = NULL;
	nline = 1 + 2*params.oversize[0] + ix_isPrim;
	MPI_Type_contiguous(nline, ntype_[0][ix_isPrim], &(ntypeSum_[0][ix_isPrim]));    //line
	MPI_Type_commit( &(ntypeSum_[0][ix_isPrim]) );
            
    }
    
} //END createType

