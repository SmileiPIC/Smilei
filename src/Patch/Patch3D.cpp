
#include "Patch3D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D constructor 
// ---------------------------------------------------------------------------------------------------------------------
Patch3D::Patch3D(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
  : Patch( params, smpi, ipatch, n_moved)
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCreation(params, smpi);
} // End Patch3D::Patch3D


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D cloning constructor 
// ---------------------------------------------------------------------------------------------------------------------
Patch3D::Patch3D(Patch3D* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
  : Patch( patch, params, smpi, ipatch, n_moved)
{
    initStep2(params);
    initStep3(params, smpi, n_moved);
    finishCloning(patch, params, smpi);
} // End Patch3D::Patch3D


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor 
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initStep2(Params& params)
{
    int xcall, ycall, zcall;
    
    // define patch coordinates
    Pcoordinates.resize(3);
    generalhilbertindexinv(params.mi[0], params.mi[1], params.mi[2], &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
#ifdef _DEBUG
    cout << "\tPatch coords : ";
    for (int iDim=0; iDim<3;iDim++)
        cout << "\t" << Pcoordinates[iDim] << " ";
    cout << endl;
#endif
    

    // 1st direction
    xcall = Pcoordinates[0]-1;
    ycall = Pcoordinates[1];
    zcall = Pcoordinates[2];
    if (params.bc_em_type_x[0]=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
    neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);
    xcall = Pcoordinates[0]+1;
    if (params.bc_em_type_x[0]=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
    neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);

    // 2st direction
    xcall = Pcoordinates[0];
    ycall = Pcoordinates[1]-1;
    zcall = Pcoordinates[2];
    if (params.bc_em_type_y[0]=="periodic" && ycall < 0) ycall += (1<<params.mi[1]);
    neighbor_[1][0] =  generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);
    ycall = Pcoordinates[1]+1;
    if (params.bc_em_type_y[0]=="periodic" && ycall >= (1<<params.mi[1])) ycall -= (1<<params.mi[1]);
    neighbor_[1][1] =  generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);

    // 3st direction
    xcall = Pcoordinates[0];
    ycall = Pcoordinates[1];
    zcall = Pcoordinates[2]-1;
    if (params.bc_em_type_z[0]=="periodic" && zcall < 0) zcall += (1<<params.mi[2]);
    neighbor_[2][0] =  generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);
    zcall = Pcoordinates[2]+1;
    if (params.bc_em_type_z[0]=="periodic" && zcall >= (1<<params.mi[2])) zcall -= (1<<params.mi[2]);
    neighbor_[2][1] =  generalhilbertindex( params.mi[0], params.mi[1], params.mi[2], xcall, ycall, zcall);

    // Corners (only used as reference during moving window operation)
    xcall = Pcoordinates[0];
    ycall = Pcoordinates[1];
    zcall = Pcoordinates[2];
#ifdef _PATCH3D_TODO
    corner_neighbor_[0][0] = 
    corner_neighbor_[0][1] = 
    corner_neighbor_[1][0] = 
    corner_neighbor_[1][1] = 
    corner_neighbor_[2][0] = 
    corner_neighbor_[2][1] = 
#endif
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI in direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initSumField( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(3);

    int patch_ndims_(3);
    int patch_nbNeighbors_(2);
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field3D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f3D->isDual_[0];
    if (field->dims_.size()>1) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f3D->isDual_[1];
        if (field->dims_.size()>2) {
            oversize2[2] *= 2;
            oversize2[2] += 1 + f3D->isDual_[2];
        }
    }

    vector<int> idx( patch_ndims_,1 );
    idx[iDim] = 0;    

    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
        std::vector<unsigned int> tmp(patch_ndims_,0);
        tmp[0] =    idx[0]  * n_elem[0] + (1-idx[0]) * oversize2[0];
        tmp[1] =    idx[1]  * n_elem[1] + (1-idx[1]) * oversize2[1];
        tmp[2] =    idx[2]  * n_elem[2] + (1-idx[2]) * oversize2[2];
        buf[iDim][iNeighbor].allocateDims( tmp );
    }
     
    int istart, ix, iy, iz;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    memset(&(idx[0]), 0, sizeof(idx[0])*idx.size());
    idx[iDim] = 1;    
        
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]][isDual[2]];
        
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = buildtag( hindex, iDim, iNeighbor );
            MPI_Isend( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, 
                       MPI_COMM_WORLD, &(f3D->MPIbuff.srequest[iDim][iNeighbor]) );
        } // END of Send
            
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1]*(buf[iDim][(iNeighbor+1)%2]).dims_[2];
            int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
            MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_3D[0][0][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, 
                       MPI_COMM_WORLD, &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
        } // END of Recv
            
    } // END for iNeighbor

} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::finalizeSumField( Field* field, int iDim )
{
    int patch_ndims_(3);
    int patch_nbNeighbors_(2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field3D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f3D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f3D->isDual_[1];
    oversize2[2] *= 2;
    oversize2[2] += 1 + f3D->isDual_[2];
    
    int istart;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
        
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &(f3D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            MPI_Wait( &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }


    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
    vector<int> idx( patch_ndims_,0 );
    idx[iDim] = 1;    
       
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
        int ix0 = idx[0]*istart;
        int iy0 = idx[1]*istart;
        int iz0 = idx[2]*istart;
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++) {
                    for (unsigned int iz=0 ; iz< (buf[iDim][(iNeighbor+1)%2]).dims_[2] ; iz++)
                        f3D->data_3D[ix0+ix][iy0+iy][iz0+iz] += (buf[iDim][(iNeighbor+1)%2])(ix,iy,iz);
                }
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
void Patch3D::initExchange( Field* field )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(3);

    int patch_ndims_(3);
    int patch_nbNeighbors_(2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);

    int istart, ix, iy, iz;

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

	vector<int> idx( patch_ndims_,0 );
	idx[iDim] = 1;

        MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]][isDual[2]];
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
                ix = idx[0]*istart;
                iy = idx[1]*istart;
                iz = idx[2]*istart;
                int tag = buildtag( hindex, iDim, iNeighbor );
                MPI_Isend( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, 
                           MPI_COMM_WORLD, &(f3D->MPIbuff.srequest[iDim][iNeighbor]) );

            } // END of Send

            if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 - (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = idx[0]*istart;
                iy = idx[1]*istart;
                iz = idx[2]*istart;
                 int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
                MPI_Irecv( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, 
                           MPI_COMM_WORLD, &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

            } // END of Recv

        } // END for iNeighbor

    } // END for iDim
} // END initExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI  (includes loop / nDim_fields_)
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::finalizeExchange( Field* field )
{
    Field3D* f3D =  static_cast<Field3D*>(field);

    int patch_ndims_(3);
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &(f3D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
             if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
               MPI_Wait( &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }

    } // END for iDim

} // END finalizeExchange( Field* field )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initExchange( Field* field, int iDim )
{
    if (field->MPIbuff.srequest.size()==0)
        field->MPIbuff.allocate(3);

    int patch_ndims_(3);
    int patch_nbNeighbors_(2);

    vector<int> idx( patch_ndims_,0 );
    idx[iDim] = 1;

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);

    int istart, ix, iy, iz;

    MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]][isDual[2]];
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

            istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = buildtag( hindex, iDim, iNeighbor );
            MPI_Isend( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, 
                       MPI_COMM_WORLD, &(f3D->MPIbuff.srequest[iDim][iNeighbor]) );

        } // END of Send

        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1- (oversize[iDim]-1) ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
            MPI_Irecv( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, 
                       MPI_COMM_WORLD, &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));

        } // END of Recv

    } // END for iNeighbor


} // END initExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::finalizeExchange( Field* field, int iDim )
{
    int patch_ndims_(3);

    Field3D* f3D =  static_cast<Field3D*>(field);

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &(f3D->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            MPI_Wait( &(f3D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
        }
    }

} // END finalizeExchange( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::createType( Params& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];
    int nz0 = params.n_space[2] + 1 + 2*params.oversize[2];
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny, nz;
    int nx_sum, ny_sum, nz_sum;
    
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;
            for (int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++) {
                nz = nz0 + iz_isPrim;
            
                // Standard Type
                ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_contiguous(params.oversize[0]*ny*nz, 
                                    MPI_DOUBLE, &(ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim]) );

                ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_vector(nx, params.oversize[1]*nz, ny, 
                                MPI_DOUBLE, &(ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim]) );

                ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_vector(nx*ny, params.oversize[2], nz, 
                                MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim]) );
            

                nx_sum = 1 + 2*params.oversize[0] + ix_isPrim;
                ny_sum = 1 + 2*params.oversize[1] + iy_isPrim;
                nz_sum = 1 + 2*params.oversize[2] + iz_isPrim;

                ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_contiguous(nx_sum*ny*nz, 
                                    MPI_DOUBLE, &(ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim]) );
            
                ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_vector(nx, ny_sum*nz, ny, 
                                MPI_DOUBLE, &(ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim]) );

                ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] = NULL;
                MPI_Type_vector(nx*ny, nz_sum, nz, 
                                MPI_DOUBLE, &(ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim]));
                MPI_Type_commit( &(ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim]) );
            
            }
        }
    }
    
} //END createType

