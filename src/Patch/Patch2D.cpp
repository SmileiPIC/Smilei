
#include "Patch2D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
    : Patch( params, smpi, domain_decomposition, ipatch, n_moved )
{
    // Test if the patch is a particle patch (Hilbert or Linearized are for VectorPatch)
    if( ( dynamic_cast<HilbertDomainDecomposition *>( domain_decomposition ) ) 
        || ( dynamic_cast<LinearizedDomainDecomposition *>( domain_decomposition ) ) ) {
        initStep2( params, domain_decomposition );
        initStep3( params, smpi, n_moved );
        finishCreation( params, smpi, domain_decomposition );
    } else { // Cartesian
        // See void Patch::set( VectorPatch& vecPatch )
        
        for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
            for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
                ntype_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntype_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntype_[2][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                
                ntype_complex_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntype_complex_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntype_complex_[2][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_complex_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_complex_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
                
            }
        }
        
    }
    
} // End Patch2D::Patch2D


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D( Patch2D *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
    : Patch( patch, params, smpi, domain_decomposition, ipatch, n_moved, with_particles )
{
    initStep2( params, domain_decomposition );
    initStep3( params, smpi, n_moved );
    finishCloning( patch, params, smpi, n_moved, with_particles );
} // End Patch2D::Patch2D


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initStep2( Params &params, DomainDecomposition *domain_decomposition )
{
    std::vector<int> xcall( 2, 0 );
    
    Pcoordinates.resize( 2 );
    Pcoordinates = domain_decomposition->getDomainCoordinates( hindex );
    
    // 1st direction
    xcall[0] = Pcoordinates[0]-1;
    xcall[1] = Pcoordinates[1];
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] < 0 ) {
        xcall[0] += domain_decomposition->ndomain_[0];
    }
    neighbor_[0][0] = domain_decomposition->getDomainId( xcall );
    
    xcall[0] = Pcoordinates[0]+1;
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] >= ( int )domain_decomposition->ndomain_[0] ) {
        xcall[0] -= domain_decomposition->ndomain_[0];
    }
    neighbor_[0][1] = domain_decomposition->getDomainId( xcall );
    
    // 2nd direction
    xcall[0] = Pcoordinates[0];
    xcall[1] = Pcoordinates[1]-1;
    if( params.EM_BCs[1][0]=="periodic" && xcall[1] < 0 ) {
        xcall[1] += domain_decomposition->ndomain_[1];
    }
    neighbor_[1][0] = domain_decomposition->getDomainId( xcall );
    
    xcall[1] = Pcoordinates[1]+1;
    if( params.EM_BCs[1][0]=="periodic" && xcall[1] >= ( int )domain_decomposition->ndomain_[1] ) {
        xcall[1] -=  domain_decomposition->ndomain_[1];
    }
    neighbor_[1][1] = domain_decomposition->getDomainId( xcall );
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ntype_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ntype_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ntype_[2][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ntypeSum_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ntypeSum_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
        }
    }
    //cout << endl;
    //cout << "Nei\t"  << "\t" << neighbor_[1][1] << endl;
    //cout << "Nei\t"  << neighbor_[0][0] << "\t" << hindex << "\t" << neighbor_[0][1] << endl;
    //cout << "Nei\t"  << "\t" << neighbor_[1][0] << endl;
    
}


Patch2D::~Patch2D()
{
    cleanType();
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI in direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initSumField( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.buf[0][0].size()==0 ) {
        field->MPIbuff.allocate( 2, field, oversize );
        
        int tagp( 0 );
        if( field->name == "Jx" ) {
            tagp = 1;
        }
        if( field->name == "Jy" ) {
            tagp = 2;
        }
        if( field->name == "Jz" ) {
            tagp = 3;
        }
        if( field->name == "Rho" ) {
            tagp = 4;
        }
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
//    int patch_ndims_(2);
    int patch_nbNeighbors_( 2 );
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D *f2D =  static_cast<Field2D *>( field );
    
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    if( field->dims_.size()>1 ) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f2D->isDual_[1];
    }
    
    int istart, ix, iy;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
    
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-iNeighbor ) * ( 0 );
            ix = ( 1-iDim )*istart;
            iy =    iDim *istart;
            int tag = f2D->MPIbuff.send_tags_[iDim][iNeighbor];
            //cout << hindex << " send to " << neighbor_[iDim][iNeighbor] << endl;
            MPI_Isend( &( ( *f2D )( ix, iy ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.srequest[iDim][iNeighbor] ) );
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            int tmp_elem = f2D->MPIbuff.buf[iDim][( iNeighbor+1 )%2].size();
            int tag = f2D->MPIbuff.recv_tags_[iDim][iNeighbor];
            //cout << hindex << " recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
            MPI_Irecv( &( f2D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
            
        } // END of Recv
        
    } // END for iNeighbor
    
} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::finalizeSumField( Field *field, int iDim )
{
    int patch_ndims_( 2 );
//    int patch_nbNeighbors_(2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D *f2D =  static_cast<Field2D *>( field );
    
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
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            //cout << hindex << " is waiting for send at " << neighbor_[iDim][iNeighbor] << endl;
            MPI_Wait( &( f2D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            //cout << hindex << " is waiting for recv from " << neighbor_[iDim][(iNeighbor+1)%2] << endl;
            MPI_Wait( &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
    int istart;
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        std::vector<unsigned int> tmp( 2, 0 );
        tmp[0] =    iDim  * n_elem[0] + ( 1-iDim ) * oversize2[0];
        if( field->dims_.size()>1 ) {
            tmp[1] = ( 1-iDim ) * n_elem[1] +    iDim  * oversize2[1];
        } else {
            tmp[1] = 1;
        }
        
        istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 );
        int ix0 = ( 1-iDim )*istart;
        int iy0 =    iDim *istart;
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            for( unsigned int ix=0 ; ix< tmp[0] ; ix++ ) {
                #pragma omp simd
                for( unsigned int iy=0 ; iy< tmp[1] ; iy++ ) {
                    ( *f2D )( ix0+ix, iy0+iy ) += f2D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][ix*tmp[1] + iy];
                }
            }
        } // END if
        
    } // END for iNeighbor
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::initExchange( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 2 );
        
        int tagp( 0 );
        if( field->name == "Bx" ) {
            tagp = 6;
        }
        if( field->name == "By" ) {
            tagp = 7;
        }
        if( field->name == "Bz" ) {
            tagp = 8;
        }
        
        if( field->name == "Ex" ) {
            tagp = 6;
        }
        if( field->name == "Ey" ) {
            tagp = 7;
        }
        if( field->name == "Ez" ) {
            tagp = 8;
        }
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
    int patch_nbNeighbors_( 2 );
    
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D *f2D =  static_cast<Field2D *>( field );
    
    int istart, ix, iy;
    
    MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = ( 1-iDim )*istart;
            iy =    iDim *istart;
            int tag = f2D->MPIbuff.send_tags_[iDim][iNeighbor];
            //cout << MPI_me_ << " Isend to " << MPI_neighbor_[iDim][iNeighbor] << " with tag " << tag << " \t name = " << field->name << endl;
            MPI_Isend( &( ( *f2D )( ix, iy ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.srequest[iDim][iNeighbor] ) );
            
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1- ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = ( 1-iDim )*istart;
            iy =    iDim *istart;
            int tag = f2D->MPIbuff.recv_tags_[iDim][iNeighbor];
            //cout << MPI_me_  << " Irecv " << MPI_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << " \t name = " << field->name << endl;
            MPI_Irecv( &( ( *f2D )( ix, iy ) ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
            
        } // END of Recv
        
    } // END for iNeighbor
    
    
} // END initExchange( Field* field, int iDim )

void Patch2D::initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 2 );
        
        int tagp( 0 );
        if( field->name == "Bx" ) {
            tagp = 6;
        }
        if( field->name == "By" ) {
            tagp = 7;
        }
        if( field->name == "Bz" ) {
            tagp = 8;
        }
        
        if( field->name == "Ex" ) {
            tagp = 6;
        }
        if( field->name == "Ey" ) {
            tagp = 7;
        }
        if( field->name == "Ez" ) {
            tagp = 8;
        }
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
    int patch_nbNeighbors_( 2 );
    
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    cField2D *f2D =  static_cast<cField2D *>( field );
    
    int istart, ix, iy;
    
    MPI_Datatype ntype = ntype_complex_[iDim][isDual[0]][isDual[1]];
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = ( 1-iDim )*istart;
            iy =    iDim *istart;
            int tag = f2D->MPIbuff.send_tags_[iDim][iNeighbor];
            //int tag = buildtag( hindex, iDim, iNeighbor, tagp );
            //cout << MPI_me_ << " Isend to " << MPI_neighbor_[iDim][iNeighbor] << " with tag " << tag << " \t name = " << field->name << endl;
            MPI_Isend( &( ( *f2D )( ix, iy ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.srequest[iDim][iNeighbor] ) );
            
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1- ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = ( 1-iDim )*istart;
            iy =    iDim *istart;
            int tag = f2D->MPIbuff.recv_tags_[iDim][iNeighbor];
            //int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, tagp );
            //cout << MPI_me_  << " Irecv " << MPI_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << " \t name = " << field->name << endl;
            MPI_Irecv( &( ( *f2D )( ix, iy ) ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
            
        } // END of Recv
        
    } // END for iNeighbor
} // END initExchangeComplex( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::finalizeExchange( Field *field, int iDim )
{
    int patch_ndims_( 2 );
    
    Field2D *f2D =  static_cast<Field2D *>( field );
    
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f2D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
} // END finalizeExchange( Field* field, int iDim )

void Patch2D::finalizeExchangeComplex( Field *field, int iDim )
{
    int patch_ndims_( 2 );
    
    cField2D *f2D =  static_cast<cField2D *>( field );
    
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f2D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f2D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
} // END finalizeExchangeComplex( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::createType( Params &params )
{
    if( ntype_[0][0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    int nx0 = params.n_space[0] + 1 + 2*oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*oversize[1];
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;
    int nx_sum, ny_sum;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            
            // Standard Type
            ntype_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( oversize[0]*ny, MPI_DOUBLE, &( ntype_[0][ix_isPrim][iy_isPrim] ) ); //line
            MPI_Type_commit( &( ntype_[0][ix_isPrim][iy_isPrim] ) );
            ntype_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_vector( nx, oversize[1], ny, MPI_DOUBLE, &( ntype_[1][ix_isPrim][iy_isPrim] ) ); // column
            MPI_Type_commit( &( ntype_[1][ix_isPrim][iy_isPrim] ) );
            
            // Still used ???
            ntype_[2][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( ny*clrw, MPI_DOUBLE, &( ntype_[2][ix_isPrim][iy_isPrim] ) ); //clrw lines
            MPI_Type_commit( &( ntype_[2][ix_isPrim][iy_isPrim] ) );
            
            ntypeSum_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            nline = 1 + 2*oversize[0] + ix_isPrim;
            //MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            
            MPI_Datatype tmpType = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( ny, MPI_DOUBLE, &( tmpType ) ); //line
            MPI_Type_commit( &( tmpType ) );
            
            MPI_Type_contiguous( nline, tmpType, &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) ); //line
            MPI_Type_commit( &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) );
            
            MPI_Type_free( &tmpType );
            
            ntypeSum_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ncol  = 1 + 2*oversize[1] + iy_isPrim;
            MPI_Type_vector( nx, ncol, ny, MPI_DOUBLE, &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) ); // column
            MPI_Type_commit( &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) );
            
            
            // Complex Type
            ntype_complex_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( 2*oversize[0]*ny, MPI_DOUBLE, &( ntype_complex_[0][ix_isPrim][iy_isPrim] ) ); //line
            MPI_Type_commit( &( ntype_complex_[0][ix_isPrim][iy_isPrim] ) );
            ntype_complex_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_vector( nx, 2*oversize[1], 2*ny, MPI_DOUBLE, &( ntype_complex_[1][ix_isPrim][iy_isPrim] ) ); // column
            MPI_Type_commit( &( ntype_complex_[1][ix_isPrim][iy_isPrim] ) );
            
            // Sum
            nx_sum = 1 + 2*oversize[0] + ix_isPrim;
            ny_sum = 1 + 2*oversize[1] + iy_isPrim;
            
            ntypeSum_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( nx_sum*ny,
                                 MPI_DOUBLE, &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) );
            MPI_Type_commit( &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) );
            
            ntypeSum_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_vector( nx, ny_sum, ny,
                             MPI_DOUBLE, &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) );
            MPI_Type_commit( &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) );
            
            // Complex sum
            ntypeSum_complex_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ntypeSum_complex_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            
        }
    }
    
} //END createType

// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::createType2( Params &params )
{
    if( ntype_[0][0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    int nx0 = params.n_space_region[0] + 1 + 2*oversize[0];
    int ny0 = params.n_space_region[1] + 1 + 2*oversize[1];
    //unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            
            // Standard Type
            ntype_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( oversize[0]*ny, MPI_DOUBLE, &( ntype_[0][ix_isPrim][iy_isPrim] ) ); //line
            MPI_Type_commit( &( ntype_[0][ix_isPrim][iy_isPrim] ) );
            ntype_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_vector( nx, oversize[1], ny, MPI_DOUBLE, &( ntype_[1][ix_isPrim][iy_isPrim] ) ); // column
            MPI_Type_commit( &( ntype_[1][ix_isPrim][iy_isPrim] ) );
            
            // Still used ???
            ntype_[2][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous(ny*params.n_space[0], MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &( ntype_[2][ix_isPrim][iy_isPrim] ) );
            
            ntypeSum_[0][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            nline = 1 + 2*oversize[0] + ix_isPrim;
            //MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            
            MPI_Datatype tmpType = MPI_DATATYPE_NULL;
            MPI_Type_contiguous( ny, MPI_DOUBLE, &( tmpType ) ); //line
            MPI_Type_commit( &( tmpType ) );
            
            MPI_Type_contiguous( nline, tmpType, &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) ); //line
            MPI_Type_commit( &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) );
            
            MPI_Type_free( &tmpType );
            
            ntypeSum_[1][ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            ncol  = 1 + 2*oversize[1] + iy_isPrim;
            MPI_Type_vector( nx, ncol, ny, MPI_DOUBLE, &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) ); // column
            MPI_Type_commit( &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) );
            
        }
    }
    
} //END createType

void Patch2D::cleanType()
{
    if( ntype_[0][0][0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            MPI_Type_free( &( ntype_[0][ix_isPrim][iy_isPrim] ) );
            MPI_Type_free( &( ntype_[1][ix_isPrim][iy_isPrim] ) );
            MPI_Type_free( &( ntype_[2][ix_isPrim][iy_isPrim] ) );
            MPI_Type_free( &( ntypeSum_[0][ix_isPrim][iy_isPrim] ) );
            MPI_Type_free( &( ntypeSum_[1][ix_isPrim][iy_isPrim] ) );
            
            if ( ntype_complex_[0][ix_isPrim][iy_isPrim] != MPI_DATATYPE_NULL ) {
                MPI_Type_free( &( ntype_complex_[0][ix_isPrim][iy_isPrim] ) );
                MPI_Type_free( &( ntype_complex_[1][ix_isPrim][iy_isPrim] ) );
                //MPI_Type_free( &(ntype_complex_[2][ix_isPrim][iy_isPrim]) );
            }
        }
    }
}


void Patch2D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
    int istart, ix, iy, iDim, iNeighbor,bufsize;
    void* b;

    bufsize = clrw*n_elem[1]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);
    iDim = 0; // We exchange only in the X direction for movewin.
    iNeighbor = 0; // We send only towards the West and receive from the East.

    MPI_Datatype ntype = ntype_[2][isDual[0]][isDual[1]]; //ntype_[2] is clrw columns.
    MPI_Status rstat    ;
    MPI_Request rrequest;


    if (MPI_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {

        istart =  2*oversize[iDim] + 1 + isDual[iDim] ;
        ix = (1-iDim)*istart;
        iy =    iDim *istart;
        MPI_Bsend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], 0, MPI_COMM_WORLD);
    } // END of Send

    //Once the message is in the buffer we can safely shift the field in memory. 
    field->shift_x(clrw);
    // and then receive the complementary field from the East.

    if (MPI_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {

        istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - clrw ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
        ix = (1-iDim)*istart;
        iy =    iDim *istart;
        MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], 0, MPI_COMM_WORLD, &rrequest);
    } // END of Recv


    if (MPI_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        MPI_Wait( &rrequest, &rstat);
    }
    MPI_Buffer_detach( &b, &bufsize);
    free(b);


} // END exchangeField_movewin

