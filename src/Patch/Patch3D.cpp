
#include "Patch3D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch3D::Patch3D( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
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
                for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                    ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    
                    ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                    ntypeSum_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                }
            }
        }
        
    }
    
} // End Patch3D::Patch3D


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch3D::Patch3D( Patch3D *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
    : Patch( patch, params, smpi, domain_decomposition, ipatch, n_moved, with_particles )
{
    initStep2( params, domain_decomposition );
    initStep3( params, smpi, n_moved );
    finishCloning( patch, params, smpi, n_moved, with_particles );
} // End Patch3D::Patch3D


// ---------------------------------------------------------------------------------------------------------------------
// Patch3D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initStep2( Params &params, DomainDecomposition *domain_decomposition )
{
    std::vector<int> xcall( 3, 0 );
    
    // define patch coordinates
    Pcoordinates.resize( 3 );
    Pcoordinates = domain_decomposition->getDomainCoordinates( hindex );
#ifdef _DEBUG
    cout << "\tPatch coords : ";
    for( int iDim=0; iDim<3; iDim++ ) {
        cout << "\t" << Pcoordinates[iDim] << " ";
    }
    cout << endl;
#endif
    
    
    // 1st direction
    xcall[0] = Pcoordinates[0]-1;
    xcall[1] = Pcoordinates[1];
    xcall[2] = Pcoordinates[2];
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] < 0 ) {
        xcall[0] += domain_decomposition->ndomain_[0];
    }
    neighbor_[0][0] = domain_decomposition->getDomainId( xcall );
    xcall[0] = Pcoordinates[0]+1;
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] >= ( int )domain_decomposition->ndomain_[0] ) {
        xcall[0] -= domain_decomposition->ndomain_[0];
    }
    neighbor_[0][1] = domain_decomposition->getDomainId( xcall );
    
    // 2st direction
    xcall[0] = Pcoordinates[0];
    xcall[1] = Pcoordinates[1]-1;
    xcall[2] = Pcoordinates[2];
    if( params.EM_BCs[1][0]=="periodic" && xcall[1] < 0 ) {
        xcall[1] += domain_decomposition->ndomain_[1];
    }
    neighbor_[1][0] =  domain_decomposition->getDomainId( xcall );
    xcall[1] = Pcoordinates[1]+1;
    if( params.EM_BCs[1][0]=="periodic" && xcall[1] >= ( int )domain_decomposition->ndomain_[1] ) {
        xcall[1] -= domain_decomposition->ndomain_[1];
    }
    neighbor_[1][1] =  domain_decomposition->getDomainId( xcall );
    
    // 3st direction
    xcall[0] = Pcoordinates[0];
    xcall[1] = Pcoordinates[1];
    xcall[2] = Pcoordinates[2]-1;
    if( params.EM_BCs[2][0]=="periodic" && xcall[2] < 0 ) {
        xcall[2] += domain_decomposition->ndomain_[2];
    }
    neighbor_[2][0] =  domain_decomposition->getDomainId( xcall );
    xcall[2] = Pcoordinates[2]+1;
    if( params.EM_BCs[2][0]=="periodic" && xcall[2] >= ( int )domain_decomposition->ndomain_[2] ) {
        xcall[2] -= domain_decomposition->ndomain_[2];
    }
    neighbor_[2][1] =  domain_decomposition->getDomainId( xcall );
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                
            }
        }
    }
    
}


Patch3D::~Patch3D()
{
    cleanType();
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI in direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initSumField( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.buf[0][0].size()==0 ) {
        field->MPIbuff.allocate( 3, field, oversize );
        
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
    
    int patch_ndims_( 3 );
    int patch_nbNeighbors_( 2 );
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D *f3D =  static_cast<Field3D *>( field );
    
    // Use a buffer per direction to exchange data before summing
    //Field3D buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f3D->isDual_[0];
    if( field->dims_.size()>1 ) {
        oversize2[1] *= 2;
        oversize2[1] += 1 + f3D->isDual_[1];
        if( field->dims_.size()>2 ) {
            oversize2[2] *= 2;
            oversize2[2] += 1 + f3D->isDual_[2];
        }
    }
    
    vector<int> idx( patch_ndims_, 1 );
    idx[iDim] = 0;
    
    int istart, ix, iy, iz;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    memset( &( idx[0] ), 0, sizeof( idx[0] )*idx.size() );
    idx[iDim] = 1;
    
    MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]][isDual[2]];
    
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-iNeighbor ) * ( 0 );
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = f3D->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( &( ( *f3D )( ix, iy, iz ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.srequest[iDim][iNeighbor] ) );
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            int tmp_elem = f3D->MPIbuff.buf[iDim][( iNeighbor+1 )%2].size();
            int tag = f3D->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &( f3D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
        } // END of Recv
        
    } // END for iNeighbor
    
} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::finalizeSumField( Field *field, int iDim )
{
    int patch_ndims_( 3 );
//    int patch_nbNeighbors_(2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D *f3D =  static_cast<Field3D *>( field );
    
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
    
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f3D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
    int istart;
    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
    vector<int> idx( 3, 1 );
    idx[iDim] = 0;
    std::vector<unsigned int> tmp( 3, 0 );
    tmp[0] =    idx[0]  * n_elem[0] + ( 1-idx[0] ) * oversize2[0];
    tmp[1] =    idx[1]  * n_elem[1] + ( 1-idx[1] ) * oversize2[1];
    tmp[2] =    idx[2]  * n_elem[2] + ( 1-idx[2] ) * oversize2[2];
    
    memset( &( idx[0] ), 0, sizeof( idx[0] )*idx.size() );
    idx[iDim] = 1;
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 );
        int ix0 = idx[0]*istart;
        int iy0 = idx[1]*istart;
        int iz0 = idx[2]*istart;
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            for( unsigned int ix=0 ; ix< tmp[0] ; ix++ ) {
                for( unsigned int iy=0 ; iy< tmp[1] ; iy++ ) {
                    for( unsigned int iz=0 ; iz< tmp[2] ; iz++ ) {
                        ( *f3D )( ix0+ix, iy0+iy, iz0+iz ) += f3D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][ix*tmp[1]*tmp[2] + iy*tmp[2] + iz];
                    }
                }
            }
        } // END if
        
    } // END for iNeighbor
    
} // END finalizeSumField


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::initExchange( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 3 );
        
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
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
    int patch_ndims_( 3 );
    int patch_nbNeighbors_( 2 );
    
    vector<int> idx( patch_ndims_, 0 );
    idx[iDim] = 1;
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D *f3D =  static_cast<Field3D *>( field );
    
    int istart, ix, iy, iz;
    
    MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]][isDual[2]];
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = f3D->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( &( ( *f3D )( ix, iy, iz ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.srequest[iDim][iNeighbor] ) );
                       
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1- ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = f3D->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &( ( *f3D )( ix, iy, iz ) ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
                       
        } // END of Recv
        
    } // END for iNeighbor
    
    
} // END initExchange( Field* field, int iDim )

void Patch3D::initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 3 );
        
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
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
    int patch_ndims_( 3 );
    int patch_nbNeighbors_( 2 );
    
    vector<int> idx( patch_ndims_, 0 );
    idx[iDim] = 1;
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    cField3D *f3D =  static_cast<cField3D *>( field );
    
    int istart, ix, iy, iz;
    
    MPI_Datatype ntype = ntype_complex_[iDim][isDual[0]][isDual[1]][isDual[2]];
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = f3D->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( &( ( *f3D )( ix, iy, iz ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.srequest[iDim][iNeighbor] ) );
                       
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1- ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = idx[0]*istart;
            iy = idx[1]*istart;
            iz = idx[2]*istart;
            int tag = f3D->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &( ( *f3D )( ix, iy, iz ) ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag,
                       MPI_COMM_WORLD, &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
                       
        } // END of Recv
        
    } // END for iNeighbor
    
    
} // END initExchangeComplex( Field* field, int iDim )

// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::finalizeExchange( Field *field, int iDim )
{
    int patch_ndims_( 3 );
    
    Field3D *f3D =  static_cast<Field3D *>( field );
    
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f3D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
} // END finalizeExchange( Field* field, int iDim )

void Patch3D::finalizeExchangeComplex( Field *field, int iDim )
{
    int patch_ndims_( 3 );
    
    cField3D *f3D =  static_cast<cField3D *>( field );
    
    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f3D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f3D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
} // END finalizeExchangeComplex( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::createType2( Params &params )
{
    if( ntype_[0][0][0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    int nx0 = params.n_space_region[0] + 1 + 2*oversize[0];
    int ny0 = params.n_space_region[1] + 1 + 2*oversize[1];
    int nz0 = params.n_space_region[2] + 1 + 2*oversize[2];
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny, nz;
    int nx_sum, ny_sum, nz_sum;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                nz = nz0 + iz_isPrim;
                
                // Standard Type
                ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_contiguous( oversize[0]*ny*nz,
                                     MPI_DOUBLE, &( ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx, oversize[1]*nz, ny*nz,
                                 MPI_DOUBLE, &( ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx*ny, oversize[2], nz,
                                 MPI_DOUBLE, &( ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );

                
                // Still used ???
                ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                if (!params.is_pxr)
                    MPI_Type_contiguous(nz*ny*params.n_space[0], MPI_DOUBLE, &(ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim]));   //clrw lines
                else
                    MPI_Type_contiguous(nz*ny*(params.n_space[0]), MPI_DOUBLE, &(ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim]));   //clrw lines
                MPI_Type_commit( &(ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim]) ); 
                
                
                nx_sum = 1 + 2*oversize[0] + ix_isPrim;
                ny_sum = 1 + 2*oversize[1] + iy_isPrim;
                nz_sum = 1 + 2*oversize[2] + iz_isPrim;
                
                ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_contiguous( nx_sum*ny*nz,
                                     MPI_DOUBLE, &( ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx, ny_sum*nz, ny*nz,
                                 MPI_DOUBLE, &( ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx*ny, nz_sum, nz,
                                 MPI_DOUBLE, &( ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
            }
        }
    }
    
}

void Patch3D::createType( Params &params )
{
    if( ntype_[0][0][0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    int nx0 = params.n_space[0] + 1 + 2*oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*oversize[1];
    int nz0 = params.n_space[2] + 1 + 2*oversize[2];
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny, nz;
    int nx_sum, ny_sum, nz_sum;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                nz = nz0 + iz_isPrim;
                
                // Standard Type
                ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_contiguous( oversize[0]*ny*nz,
                                     MPI_DOUBLE, &( ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx, oversize[1]*nz, ny*nz,
                                 MPI_DOUBLE, &( ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx*ny, oversize[2], nz,
                                 MPI_DOUBLE, &( ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                //// Complex Type V0
                //ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                //MPI_Type_contiguous(params.oversize[0]*ny*nz,
                //                    MPI_CXX_DOUBLE_COMPLEX, &(ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim]));
                //MPI_Type_commit( &(ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim]) );
                //
                //ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                //MPI_Type_vector(nx, params.oversize[1]*nz, ny*nz,
                //                MPI_CXX_DOUBLE_COMPLEX, &(ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim]));
                //MPI_Type_commit( &(ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim]) );
                //
                //ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                //MPI_Type_vector(nx*ny, params.oversize[2], nz,
                //                MPI_CXX_DOUBLE_COMPLEX, &(ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim]));
                //MPI_Type_commit( &(ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim]) );
                
                // Complex Type
                ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_contiguous( 2*oversize[0]*ny*nz,
                                     MPI_DOUBLE, &( ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx, oversize[1]*2*nz, ny*2*nz,
                                 MPI_DOUBLE, &( ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx*ny, 2*oversize[2], 2*nz,
                                 MPI_DOUBLE, &( ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                nx_sum = 1 + 2*oversize[0] + ix_isPrim;
                ny_sum = 1 + 2*oversize[1] + iy_isPrim;
                nz_sum = 1 + 2*oversize[2] + iz_isPrim;
                
                // Standard sum
                ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_contiguous( nx_sum*ny*nz,
                                     MPI_DOUBLE, &( ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx, ny_sum*nz, ny*nz,
                                 MPI_DOUBLE, &( ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                MPI_Type_vector( nx*ny, nz_sum, nz,
                                 MPI_DOUBLE, &( ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_commit( &( ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                // Complex sum
                ntypeSum_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                ntypeSum_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                
            }
        }
    }
    
} //END createType

void Patch3D::cleanType()
{
    if( ntype_[0][0][0][0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                MPI_Type_free( &( ntype_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_free( &( ntype_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_free( &( ntype_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                if (ntype_[3][0][0][0] != MPI_DATATYPE_NULL) 
                    MPI_Type_free( &(ntype_[3][ix_isPrim][iy_isPrim][iz_isPrim]) );
                MPI_Type_free( &( ntypeSum_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_free( &( ntypeSum_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                MPI_Type_free( &( ntypeSum_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                
                if ( ntype_complex_[0][0][0][0] != MPI_DATATYPE_NULL ) {
                    MPI_Type_free( &( ntype_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                    MPI_Type_free( &( ntype_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                    MPI_Type_free( &( ntype_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim] ) );
                    //MPI_Type_free( &(ntypeSum_complex_[0][ix_isPrim][iy_isPrim][iz_isPrim]) );
                    //MPI_Type_free( &(ntypeSum_complex_[1][ix_isPrim][iy_isPrim][iz_isPrim]) );
                    //MPI_Type_free( &(ntypeSum_complex_[2][ix_isPrim][iy_isPrim][iz_isPrim]) );
                }
            }
        }
    }
}


void Patch3D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field3D* f3D =  static_cast<Field3D*>(field);
    int istart, ix, iy, iz, iDim, iNeighbor,bufsize;
    void* b;

    bool pxr(false);
    if ( (!isDual[0])&&(!isDual[1])&&(!isDual[2]) )
        pxr = true;

    if (!pxr)
        bufsize = clrw*n_elem[1]*n_elem[2]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    else
        bufsize = (clrw+oversize[0]+1)*n_elem[1]*n_elem[2]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);
    iDim = 0; // We exchange only in the X direction for movewin.
    iNeighbor = 0; // We send only towards the West and receive from the East.

    MPI_Datatype ntype = ntype_[3][isDual[0]][isDual[1]][isDual[2]]; //ntype_[3] is clrw plans.
    MPI_Status rstat    ;
    MPI_Request rrequest;


    int patch_ndims_(3);
    vector<int> idx( patch_ndims_,0 );
    idx[iDim] = 1;

    if (MPI_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
        if (!pxr)
            istart =  2*oversize[iDim] + 1 + isDual[iDim] ;
        else
            istart =  oversize[iDim];
        ix = idx[0]*istart;
        iy = idx[1]*istart;
        iz = idx[2]*istart;
        MPI_Bsend( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], 0, MPI_COMM_WORLD);
    } // END of Send

    //Once the message is in the buffer we can safely shift the field in memory. 
    field->shift_x(clrw);
    // and then receive the complementary field from the East.

    if (MPI_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
         if (!pxr)
             istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - clrw ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
         else
             istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - clrw - oversize[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
        ix = idx[0]*istart;
        iy = idx[1]*istart;
        iz = idx[2]*istart;
        MPI_Irecv( &(f3D->data_3D[ix][iy][iz]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], 0, MPI_COMM_WORLD, &rrequest);
    } // END of Recv


    if (MPI_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        MPI_Wait( &rrequest, &rstat);
    }
    MPI_Buffer_detach( &b, &bufsize);
    free(b);


} // END exchangeField_movewin
