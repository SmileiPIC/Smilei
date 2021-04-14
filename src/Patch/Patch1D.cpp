
#include "Patch1D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
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
            ntype_[ix_isPrim] = MPI_DATATYPE_NULL;
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
        ntype_[ix_isPrim] = MPI_DATATYPE_NULL;
    }
    
}

Patch1D::~Patch1D()
{
    cleanType();
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch1D::createType2( Params &params )
{
    if( ntype_[0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    unsigned int clrw = params.clrw;
    
    // int ny = oversize[0];
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
    
        ntype_[ix_isPrim] = MPI_DATATYPE_NULL;
        MPI_Type_contiguous( clrw, MPI_DOUBLE, &( ntype_[ix_isPrim] ) ); //clrw lines
        MPI_Type_commit( &( ntype_[ix_isPrim] ) );
        
    }
}

void Patch1D::cleanType()
{
    if( ntype_[0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        MPI_Type_free( &( ntype_[ix_isPrim] ) );
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
