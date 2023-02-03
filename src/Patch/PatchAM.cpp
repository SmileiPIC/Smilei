
#include "PatchAM.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "Hilbert_functions.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// PatchAM constructor
// ---------------------------------------------------------------------------------------------------------------------
PatchAM::PatchAM( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
    : Patch( params, smpi, domain_decomposition, ipatch )
{
    // Test if the patch is a particle patch (Hilbert or Linearized are for VectorPatch)
    if( ( dynamic_cast<HilbertDomainDecomposition *>( domain_decomposition ) )
        || ( dynamic_cast<LinearizedDomainDecomposition *>( domain_decomposition ) ) ) {
        initStep2( params, domain_decomposition );
        initInvR( params );
        initStep3( params, smpi, n_moved );
        finishCreation( params, smpi, domain_decomposition );
    } else { // Cartesian
        // See void Patch::set( VectorPatch& vecPatch )
        
        for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
            for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
                ntype_complex_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            }
        }
        
    }
    
} // End PatchAM::PatchAM


// ---------------------------------------------------------------------------------------------------------------------
// PatchAM cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
PatchAM::PatchAM( PatchAM *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
    : Patch( patch, params, smpi, ipatch )
{
    initStep2( params, domain_decomposition );
    initInvR( params );
    initStep3( params, smpi, n_moved );
    finishCloning( patch, params, smpi, n_moved, with_particles );
} // End PatchAM::PatchAM


// ---------------------------------------------------------------------------------------------------------------------
// PatchAM second initializer :
//   - Pcoordinates, neighbor_ resized in Patch constructor
// ---------------------------------------------------------------------------------------------------------------------
void PatchAM::initStep2( Params &params, DomainDecomposition *domain_decomposition )
{
    Pcoordinates = domain_decomposition->getDomainCoordinates( hindex );
    
    std::vector<int> xcall( 2, 0 );
    
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
            ntype_complex_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            
        }
    }
    
}

void PatchAM::initInvR( Params &params ) {
    int j_glob_ = Pcoordinates[1]*size_[1]-oversize[1]; //cell_starting_global_index is only define later during patch creation.
    int nr_p = size_[1]+1+2*oversize[1];
    double dr = params.cell_length[1];
    invR.resize( nr_p );

    if (!params.is_spectral){
        invRd.resize( nr_p+1 );
        for( int j = 0; j< nr_p; j++ ) {
            if( j_glob_ + j == 0 ) {
                invR[j] = 8./dr; // No Verboncoeur correction
                //invR[j] = 64./(13.*dr); // Order 2 Verboncoeur correction
            } else {
                invR[j] = 1./abs(((double)j_glob_ + (double)j)*dr);
            }
        }
        for( int j = 0; j< nr_p + 1; j++ ) {
            invRd[j] = 1./abs(((double)j_glob_ + (double)j - 0.5)*dr);
        }
     } else { // if spectral, primal grid shifted by half cell length
        for( int j = 0; j< nr_p; j++ ) {
            invR[j] = 1./( ( (double)j_glob_ +(double)j + 0.5)*dr);
        }
     }
}

PatchAM::~PatchAM()
{
    cleanType();
}

// ---------------------------------------------------------------------------------------------------------------------
// Initialize current patch sum Fields communications through MPI in direction iDim
// Intra-MPI process communications managed by memcpy in SyncVectorPatch::sum()
// ---------------------------------------------------------------------------------------------------------------------
void PatchAM::initSumFieldComplex( Field *field, int iDim, SmileiMPI *smpi )
{
    if( field->MPIbuff.srequest.size()==0 ) {
        field->MPIbuff.allocate( nDim_fields_ );
        
        int tagp( 0 );
        if( field->name == "Jl" ) {
            tagp = 1;
        }
        if( field->name == "Jr" ) {
            tagp = 2;
        }
        if( field->name == "Jt" ) {
            tagp = 3;
        }
        if( field->name == "Rho" ) {
            tagp = 4;
        }
        
        field->MPIbuff.defineTags( this, smpi, tagp );
    }
    
    int patch_nbNeighbors_( 2 );
    
    for( int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++ ) {
    
        if( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            int tag = field->MPIbuff.send_tags_[iDim][iNeighbor];
            MPI_Isend( static_cast<cField*>(field->sendFields_[iDim*2+iNeighbor])->cdata_, 2*field->sendFields_[iDim*2+iNeighbor]->number_of_points_, MPI_DOUBLE, MPI_neighbor_[iDim][iNeighbor], tag,
                       MPI_COMM_WORLD, &( field->MPIbuff.srequest[iDim][iNeighbor] ) );
        } // END of Send
        
        if( is_a_MPI_neighbor( iDim, ( iNeighbor+1 )%2 ) ) {
            int tag = field->MPIbuff.recv_tags_[iDim][iNeighbor];
            MPI_Irecv( static_cast<cField*>(field->recvFields_[iDim*2+(iNeighbor+1)%2])->cdata_, 2*field->recvFields_[iDim*2+(iNeighbor+1)%2]->number_of_points_, MPI_DOUBLE, MPI_neighbor_[iDim][( iNeighbor+1 )%2], tag,
                       MPI_COMM_WORLD, &( field->MPIbuff.rrequest[iDim][( iNeighbor+1 )%2] ) );

        } // END of Recv
        
    } // END for iNeighbor
    
} // END initSumFieldComplex

// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void PatchAM::createType2( Params &params )
{
    if( ntype_complex_[0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    // int nx0 = params.region_size_[0] + 1 + 2*oversize[0];
    int ny0 = params.region_size_[1] + 1 + 2*oversize[1];
    //unsigned int clrw = params.cluster_width_;

    // MPI_Datatype ntype_[nDim][primDual][primDual]
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        // int nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            int ny = ny0 + iy_isPrim;
            
            // Still used ??? Yes, for moving window and SDMD
            ntype_complex_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous(2*ny*params.patch_size_[0], MPI_DOUBLE, &(ntype_complex_[ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &( ntype_complex_[ix_isPrim][iy_isPrim] ) );
        }
    }
    
} //END createType

void PatchAM::cleanType()
{
}

void PatchAM::exchangeField_movewin( Field* field, int nshift )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    cField2D* f2D =  static_cast<cField2D*>(field);

    int bufsize = 2*nshift*n_elem[1]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    void* b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);

    MPI_Status rstat    ;
    MPI_Request rrequest;

    if( MPI_neighbor_[0][0]!=MPI_PROC_NULL ) {
        int ix = 2*oversize[0] + 1 + isDual[0];
        int iy =   0;
        MPI_Bsend(  &( ( *f2D )( ix, iy ) ), 2*nshift*n_elem[1], MPI_DOUBLE, MPI_neighbor_[0][0], 0, MPI_COMM_WORLD);
    } // END of Send

    //Once the message is in the buffer we can safely shift the field in memory.
    field->shift_x(nshift);
    // and then receive the complementary field from the East.

    if( MPI_neighbor_[0][1]!=MPI_PROC_NULL ) {
        int ix = n_elem[0] - nshift;
        int iy =   0 ;
        MPI_Irecv(  &( ( *f2D )( ix, iy ) ), 2*nshift*n_elem[1], MPI_DOUBLE, MPI_neighbor_[0][1], 0, MPI_COMM_WORLD, &rrequest);
        MPI_Wait( &rrequest, &rstat);
    }
    MPI_Buffer_detach( &b, &bufsize);
    free(b);

} // END exchangeField_movewin


void PatchAM::computePoynting() {
    if( isBoundary( 0, 0 ) ) {
        EMfields->computePoynting( 0, 0 );
    }
    if( isBoundary( 0, 1 ) ) {
        EMfields->computePoynting( 0, 1 );
    }
    if( isBoundary( 1, 1 ) ) {
        EMfields->computePoynting( 1, 1 );
    }
    // No poynting from axis
    EMfields->poynting_inst[0][1] = 0;
    EMfields->poynting[0][1] = 0;
}
