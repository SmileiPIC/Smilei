
#include "Patch2D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D constructor 
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
  : Patch( params, smpi, ipatch, n_moved)
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCreation(params, smpi);
} // End Patch2D::Patch2D


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D cloning constructor 
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D(Patch2D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved, bool with_particles = true)
    : Patch( patch, params, smpi, ipatch, n_moved, with_particles )
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCloning(patch, params, smpi, with_particles);
} // End Patch2D::Patch2D


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor 
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initStep2(Params& params)
{
    int xcall, ycall;
    
    Pcoordinates.resize(2);
    generalhilbertindexinv(params.mi[0], params.mi[1], &Pcoordinates[0], &Pcoordinates[1], hindex);
    
    // 1st direction
    xcall = Pcoordinates[0]-1;
    ycall = Pcoordinates[1];
    if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
    neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    xcall = Pcoordinates[0]+1;
    if (params.bc_em_type_x[0]=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
    neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    
    // 2nd direction
    xcall = Pcoordinates[0];
    ycall = Pcoordinates[1]-1;
    if (params.bc_em_type_y[0]=="periodic" && ycall < 0) ycall += (1<<params.mi[1]);
    neighbor_[1][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    ycall = Pcoordinates[1]+1;
    if (params.bc_em_type_y[0]=="periodic" && ycall >= (1<<params.mi[1])) ycall -= (1<<params.mi[1]);
    neighbor_[1][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    
    // Corners
    //xcall = Pcoordinates[0]+1;
    //if (params.bc_em_type_x[0]=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
    //corner_neighbor_[1][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    //xcall = Pcoordinates[0]-1;
    //if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
    //corner_neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    //ycall = Pcoordinates[1]-1;
    //if (params.bc_em_type_y[0]=="periodic" && ycall < 0) ycall += (1<<params.mi[1]);
    //corner_neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
    //xcall = Pcoordinates[0]+1;
    //if (params.bc_em_type_x[0]=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
    //corner_neighbor_[1][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
}


Patch2D::~Patch2D()
{
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            MPI_Type_free( &(ntype_[0][ix_isPrim][iy_isPrim]) );
            MPI_Type_free( &(ntype_[1][ix_isPrim][iy_isPrim]) );
            MPI_Type_free( &(ntype_[2][ix_isPrim][iy_isPrim]) );
            MPI_Type_free( &(ntypeSum_[0][ix_isPrim][iy_isPrim]) );
            MPI_Type_free( &(ntypeSum_[1][ix_isPrim][iy_isPrim]) );            
        }
    }
}


void Patch2D::reallyinitSumField( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(2);

    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    if (field->dims_.size()>1) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f2D->isDual_[1];
    }
    
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
        std::vector<unsigned int> tmp(patch_ndims_,0);
        tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
        if (field->dims_.size()>1)
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
        else 
            tmp[1] = 1;
        buf[iDim][iNeighbor].allocateDims( tmp );
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI in direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initSumField( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(2);

//    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    if (field->dims_.size()>1) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f2D->isDual_[1];
    }
    
    /*for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
        std::vector<unsigned int> tmp(patch_ndims_,0);
        tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
        if (field->dims_.size()>1)
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
        else 
            tmp[1] = 1;
        buf[iDim][iNeighbor].allocateDims( tmp );
    }*/
     
    int istart, ix, iy;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
        
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
        
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
            ix = (1-iDim)*istart;
            iy =    iDim *istart;
            int tag = send_tags_[iDim][iNeighbor];
            //cout << hindex << " send to " << neighbor_[iDim][iNeighbor] << endl;
            //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );
            MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );
        } // END of Send
            
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
            int tag = recv_tags_[iDim][iNeighbor];
            //cout << hindex << " recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
            //MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
            MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
        } // END of Recv
            
    } // END for iNeighbor

} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::finalizeSumField( Field* field, int iDim )
{
    int patch_ndims_(2);
//    int patch_nbNeighbors_(2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f2D->isDual_[1];
    
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
        
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            //cout << hindex << " is waiting for send at " << neighbor_[iDim][iNeighbor] << endl;
            MPI_Wait( &(f2D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            //cout << hindex << " is waiting for recv from " << neighbor_[iDim][(iNeighbor+1)%2] << endl;        
            MPI_Wait( &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }
}

void Patch2D::reallyfinalizeSumField( Field* field, int iDim )
{
//    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f2D->isDual_[1];

    int istart;
    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
       
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
        int ix0 = (1-iDim)*istart;
        int iy0 =    iDim *istart;
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                    f2D->data_2D[ix0+ix][iy0+iy] += (buf[iDim][(iNeighbor+1)%2])(ix,iy);
            }
        } // END if
            
    } // END for iNeighbor
        

    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
        buf[iDim][iNeighbor].deallocateDims();
    }

} // END finalizeSumField


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initExchange( Field* field )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(2);

    int patch_ndims_(2);
    int patch_nbNeighbors_(2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    int istart, ix, iy;

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

        MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                int tag = send_tags_[iDim][iNeighbor];
                //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );
                MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );

            } // END of Send

            if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 - (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                int tag = recv_tags_[iDim][iNeighbor];
                //MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));
                MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

            } // END of Recv

        } // END for iNeighbor

    } // END for iDim
} // END initExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI  (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::finalizeExchange( Field* field )
{
    Field2D* f2D =  static_cast<Field2D*>(field);

    int patch_ndims_(2);
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &(f2D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
             if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
               MPI_Wait( &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }

    } // END for iDim

} // END finalizeExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initExchange( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(2);

    int patch_nbNeighbors_(2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    int istart, ix, iy;

    MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

            istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = (1-iDim)*istart;
            iy =    iDim *istart;
            int tag = send_tags_[iDim][iNeighbor];
            //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );
            MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.srequest[iDim][iNeighbor]) );

        } // END of Send

        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1- (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
            ix = (1-iDim)*istart;
            iy =    iDim *istart;
            int tag = recv_tags_[iDim][iNeighbor];
            //MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));
            MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

        } // END of Recv

    } // END for iNeighbor


} // END initExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::finalizeExchange( Field* field, int iDim )
{
    int patch_ndims_(2);

    Field2D* f2D =  static_cast<Field2D*>(field);

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &(f2D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            MPI_Wait( &(f2D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }

} // END finalizeExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::createType( Params& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;
    
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;
            
            // Standard Type
            ntype_[0][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(params.oversize[0]*ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntype_[0][ix_isPrim][iy_isPrim]) );
            ntype_[1][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_vector(nx, params.oversize[1], ny, MPI_DOUBLE, &(ntype_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntype_[1][ix_isPrim][iy_isPrim]) );
            
            // Still used ???
            ntype_[2][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(ny*clrw, MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &(ntype_[2][ix_isPrim][iy_isPrim]) );
            
            ntypeSum_[0][ix_isPrim][iy_isPrim] = NULL;
            nline = 1 + 2*params.oversize[0] + ix_isPrim;
            //MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            
            MPI_Datatype tmpType = NULL;
            MPI_Type_contiguous(ny, MPI_DOUBLE, &(tmpType));    //line
            MPI_Type_commit( &(tmpType) );
            
            MPI_Type_contiguous(nline, tmpType, &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntypeSum_[0][ix_isPrim][iy_isPrim]) );
            
            MPI_Type_free( &tmpType );
            
            ntypeSum_[1][ix_isPrim][iy_isPrim] = NULL;
            ncol  = 1 + 2*params.oversize[1] + iy_isPrim;
            MPI_Type_vector(nx, ncol, ny, MPI_DOUBLE, &(ntypeSum_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntypeSum_[1][ix_isPrim][iy_isPrim]) );
            
        }
    }
    
} //END createType

