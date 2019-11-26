
#include "Patch1D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch1D constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch1D::Patch1D( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
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
            ntype_[0][ix_isPrim] = MPI_DATATYPE_NULL;
            ntype_[1][ix_isPrim] = MPI_DATATYPE_NULL;
            ntypeSum_[0][ix_isPrim] = MPI_DATATYPE_NULL;
            
            ntype_complex_[0][ix_isPrim] = MPI_DATATYPE_NULL;
            ntype_complex_[1][ix_isPrim] = MPI_DATATYPE_NULL;
            
            ntypeSum_complex_[0][ix_isPrim] = MPI_DATATYPE_NULL;
            
        }
        
    }
    
} // End Patch1D::Patch1D

// ---------------------------------------------------------------------------------------------------------------------
// Patch1D cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch1D::Patch1D( Patch1D *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
    : Patch( patch, params, smpi, domain_decomposition, ipatch, n_moved, with_particles )
{
    initStep2( params, domain_decomposition );
    initStep3( params, smpi, n_moved );
    finishCloning( patch, params, smpi, n_moved, with_particles );
} // End Patch1D::Patch1D


// ---------------------------------------------------------------------------------------------------------------------
// Patch1D second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initStep2( Params &params, DomainDecomposition *domain_decomposition )
{
    std::vector<int> xcall( 1, 0 );
    
    Pcoordinates[0] = hindex;
    
    // 1st direction
    xcall[0] = Pcoordinates[0]-1;
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] < 0 ) {
        xcall[0] += domain_decomposition->ndomain_[0];
    }
    neighbor_[0][0] = domain_decomposition->getDomainId( xcall );
    xcall[0] = Pcoordinates[0]+1;
    if( params.EM_BCs[0][0]=="periodic" && xcall[0] >= ( int )domain_decomposition->ndomain_[0] ) {
        xcall[0] -= domain_decomposition->ndomain_[0];
    }
    neighbor_[0][1] = domain_decomposition->getDomainId( xcall );
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        ntype_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        ntype_[1][ix_isPrim] = MPI_DATATYPE_NULL;
        ntypeSum_[0][ix_isPrim] = MPI_DATATYPE_NULL;
    }
    
}

Patch1D::~Patch1D()
{
    cleanType();
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initSumField( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.buf[0][0].size()==0 ) {
        field->MPIbuff.allocate( 1, field, oversize );
        
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
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D *f1D =  static_cast<Field1D *>( field );
    
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
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-iNeighbor ) * ( 0 );
            ix = ( 1-iDim )*istart;
            int tag = f1D->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( &( f1D->data_[ix] ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.srequest[iDim][iNeighbor] ) );
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            int tmp_elem = f1D->MPIbuff.buf[iDim][( iNeighbor+1 )%2].size();
            int tag = f1D->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &( f1D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
        } // END of Recv
        
    } // END for iNeighbor
    
} // END initSumField


// ---------------------------------------------------------------------------------------------------------------------
// Finalize current patch sum Fields communications through MPI for direction iDim
// Proceed to the local reduction
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::finalizeSumField( Field *field, int iDim )
{
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D *f1D =  static_cast<Field1D *>( field );
    
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
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f1D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        int istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 );
        int ix0 = ( 1-iDim )*istart;
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            unsigned int tmp = f1D->MPIbuff.buf[iDim][( iNeighbor+1 )%2].size();
            for( unsigned int ix=0 ; ix<tmp ; ix++ ) {
                f1D->data_[ix0+ix] += f1D->MPIbuff.buf[iDim][( iNeighbor+1 )%2][ix];
            }
        } // END if
    } // END for iNeighbor
}


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::initExchange( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 1 );
        
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
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D *f1D =  static_cast<Field1D *>( field );
    
    int istart, ix;
    
    MPI_Datatype ntype = ntype_[iDim][isDual[0]];
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = ( 1-iDim )*istart;
            int tag = f1D->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( &( f1D->data_[ix] ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.srequest[iDim][iNeighbor] ) );
            
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1 - ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = ( 1-iDim )*istart;
            int tag = f1D->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( &( f1D->data_[ix] ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
            
        } // END of Recv
        
    } // END for iNeighbor
    
} // END initExchange( Field* field, int iDim )

void Patch1D::initExchangeComplex( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( 1 );
        
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
    
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    cField1D *f1D =  static_cast<cField1D *>( field );
    
    int istart, ix;
    
    MPI_Datatype ntype = ntype_complex_[iDim][isDual[0]];
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
        
            istart = iNeighbor * ( n_elem[iDim]- ( 2*oversize[iDim]+1+isDual[iDim] ) ) + ( 1-iNeighbor ) * ( oversize[iDim] + 1 + isDual[iDim] );
            ix = ( 1-iDim )*istart;
            int tag = f1D->MPIbuff.send_tags_[iDim][iNeighbor];
            //MPI_Isend( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.srequest[iDim][iNeighbor]) );
            MPI_Isend( &( ( *f1D )( ix ) ), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.srequest[iDim][iNeighbor] ) );
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
        
            istart = ( ( iNeighbor+1 )%2 ) * ( n_elem[iDim] - 1 - ( oversize[iDim]-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
            ix = ( 1-iDim )*istart;
            int tag = f1D->MPIbuff.recv_tags_[iDim][iNeighbor];
            //MPI_Irecv( &(f1D->data_[ix]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f1D->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]));
            MPI_Irecv( &( ( *f1D )( ix ) ), 1, ntype, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag, MPI_COMM_WORLD, &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );
        } // END of Recv
        
    } // END for iNeighbor
} // END initExchangeComplex( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch exhange Fields communications through MPI for direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::finalizeExchange( Field *field, int iDim )
{
    Field1D *f1D =  static_cast<Field1D *>( field );
    
    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f1D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
} // END finalizeExchange( Field* field, int iDim )

void Patch1D::finalizeExchangeComplex( Field *field, int iDim )
{
    cField1D *f1D =  static_cast<cField1D *>( field );
    
    MPI_Status sstat    [nDim_fields_][2];
    MPI_Status rstat    [nDim_fields_][2];
    
    for( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            MPI_Wait( &( f1D->MPIbuff.srequest[iDim][iNeighbor] ), &( sstat[iDim][iNeighbor] ) );
        }
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            MPI_Wait( &( f1D->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ), &( rstat[iDim][( iNeighbor+1 )%2] ) );
        }
    }
    
} // END finalizeExchangeComplex( Field* field, int iDim )


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::createType2( Params &params )
{
    if( ntype_[0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual]
    int ny = oversize[0];
    int nline;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
    
        // Standard Type
        ntype_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( ny, MPI_DOUBLE, &( ntype_[0][ix_isPrim] ) ); //line
        MPI_Type_commit( &( ntype_[0][ix_isPrim] ) );
        
        ntype_[1][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( clrw, MPI_DOUBLE, &( ntype_[1][ix_isPrim] ) ); //clrw lines
        MPI_Type_commit( &( ntype_[1][ix_isPrim] ) );
        
        ntypeSum_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        
        MPI_Datatype tmpType = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( 1, MPI_DOUBLE, &( tmpType ) ); //line
        MPI_Type_commit( &( tmpType ) );
        
        
        nline = 1 + 2*params.oversize[0] + ix_isPrim;
        MPI_Type_contiguous( nline, tmpType, &( ntypeSum_[0][ix_isPrim] ) ); //line
        MPI_Type_commit( &( ntypeSum_[0][ix_isPrim] ) );
        
        MPI_Type_free( &tmpType );
        
    }
}

void Patch1D::createType( Params &params )
{
    if( ntype_[0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual]
    int ny = oversize[0];
    int nline;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
    
        // Standard Type
        ntype_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( ny, MPI_DOUBLE, &( ntype_[0][ix_isPrim] ) ); //line
        MPI_Type_commit( &( ntype_[0][ix_isPrim] ) );
        
        ntype_[1][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( clrw, MPI_DOUBLE, &( ntype_[1][ix_isPrim] ) ); //clrw lines
        MPI_Type_commit( &( ntype_[1][ix_isPrim] ) );
        
        ntypeSum_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        
        MPI_Datatype tmpType = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( 1, MPI_DOUBLE, &( tmpType ) ); //line
        MPI_Type_commit( &( tmpType ) );
        
        
        nline = 1 + 2*params.oversize[0] + ix_isPrim;
        MPI_Type_contiguous( nline, tmpType, &( ntypeSum_[0][ix_isPrim] ) ); //line
        MPI_Type_commit( &( ntypeSum_[0][ix_isPrim] ) );
        
        MPI_Type_free( &tmpType );
        
        
        // Complex type
        ntype_complex_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( 2*ny, MPI_DOUBLE, &( ntype_complex_[0][ix_isPrim] ) ); //line
        MPI_Type_commit( &( ntype_complex_[0][ix_isPrim] ) );
        
        ntype_complex_[1][ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( 2*clrw, MPI_DOUBLE, &( ntype_complex_[1][ix_isPrim] ) ); //clrw lines
        MPI_Type_commit( &( ntype_complex_[1][ix_isPrim] ) );
        
        ntypeSum_complex_[0][ix_isPrim] = MPI_DATATYPE_NULL;
        
        // MPI_Datatype tmpTypeComplex = MPI_DATATYPE_NULL;
        // MPI_Type_contiguous(2, MPI_DOUBLE, &(tmpTypeComplex));    //line
        // MPI_Type_commit( &(tmpTypeComplex) );
        //
        //
        // nline = 1 + 2*params.oversize[0] + ix_isPrim;
        // MPI_Type_contiguous(nline, tmpType, &(ntypeSum_[0][ix_isPrim]));    //line
        // MPI_Type_commit( &(ntypeSum_[0][ix_isPrim]) );
        //
        // MPI_Type_free( &tmpTypeComplex );
        
        
        
        
        
    }
    
} //END createType

void Patch1D::cleanType()
{
    if( ntype_[0][0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        MPI_Type_free( &( ntype_[0][ix_isPrim] ) );
        MPI_Type_free( &( ntype_[1][ix_isPrim] ) );
        MPI_Type_free( &( ntypeSum_[0][ix_isPrim] ) );

        if ( ntype_complex_[0][ix_isPrim] != MPI_DATATYPE_NULL ) {
            MPI_Type_free( &( ntype_complex_[0][ix_isPrim] ) );
            MPI_Type_free( &( ntype_complex_[1][ix_isPrim] ) );
            //MPI_Type_free( &(ntypeSum_complex_[0][ix_isPrim]) );
        }
    }
}

void Patch1D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
    int istart, /*iDim,*/ iNeighbor, bufsize;
    void* b;

//    iDim = 0; // We exchange only in the X direction for movewin.
    iNeighbor = 0; // We send only towards the West and receive from the East.

    bufsize = clrw * sizeof(double) + 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);

    // Loop over dimField
    // See sumField for details
    MPI_Status rstat;
    MPI_Request rrequest;


    if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
        istart = 2*oversize[0] + 1 + isDual[0]  ;
        MPI_Bsend( &(f1D->data_[istart]), clrw, MPI_DOUBLE, neighbor_[0][iNeighbor], 0, MPI_COMM_WORLD);
    } // END of Send

    field->shift_x(clrw);

    if (MPI_neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[0] - clrw ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
        MPI_Irecv( &(f1D->data_[istart]), clrw, MPI_DOUBLE, MPI_neighbor_[0][(iNeighbor+1)%2], 0, MPI_COMM_WORLD, &rrequest );
    } // END of Recv


    if (MPI_neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        MPI_Wait( &rrequest, &rstat );
    }

    MPI_Buffer_detach( &b, &bufsize);
    free(b);


} // END exchangeField
