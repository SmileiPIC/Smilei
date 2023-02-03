
#include "Patch2D.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "DomainDecompositionFactory.h"
#include "Species.h"
#include "Particles.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved )
    : Patch( params, smpi, domain_decomposition, ipatch )
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
                ntype_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            }
        }
        
    }
    
} // End Patch2D::Patch2D


// ---------------------------------------------------------------------------------------------------------------------
// Patch2D cloning constructor
// ---------------------------------------------------------------------------------------------------------------------
Patch2D::Patch2D( Patch2D *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true )
    : Patch( patch, params, smpi, ipatch )
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
            ntype_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
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
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch2D::createType2( Params &params )
{
    if( ntype_[0][0] != MPI_DATATYPE_NULL ) {
        return;
    }
    
    //sint nx0 = params.region_size_[0] + 1 + 2*oversize[0];
    int ny0 = params.region_size_[1] + 1 + 2*oversize[1];
    //unsigned int clrw = params.cluster_width_;
    
    int ny;
    //int nline, ncol;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        //nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            
            // Still used ???
            ntype_[ix_isPrim][iy_isPrim] = MPI_DATATYPE_NULL;
            MPI_Type_contiguous(ny*params.patch_size_[0], MPI_DOUBLE, &(ntype_[ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &( ntype_[ix_isPrim][iy_isPrim] ) );

        }
    }
    
} //END createType

void Patch2D::cleanType()
{
    if( ntype_[0][0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            MPI_Type_free( &( ntype_[ix_isPrim][iy_isPrim] ) );
        }
    }
}


void Patch2D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    int bufsize = clrw*n_elem[1]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    void *b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);

    MPI_Status rstat    ;
    MPI_Request rrequest;

    if( MPI_neighbor_[0][0]!=MPI_PROC_NULL ) {
        int ix = 2*oversize[0] + 1 + isDual[0];
        int iy =  0 ;
        MPI_Bsend( &(f2D->data_2D[ix][iy]), clrw*n_elem[1], MPI_DOUBLE, MPI_neighbor_[0][0], 0, MPI_COMM_WORLD);
    } // END of Send

    //Once the message is in the buffer we can safely shift the field in memory.
    field->shift_x(clrw);
    // and then receive the complementary field from the East.

    if( MPI_neighbor_[0][1]!=MPI_PROC_NULL ) {
        int ix = n_elem[0] - clrw  ;
        int iy =  0 ;
        MPI_Irecv( &(f2D->data_2D[ix][iy]), clrw*n_elem[1], MPI_DOUBLE, MPI_neighbor_[0][1], 0, MPI_COMM_WORLD, &rrequest);
        MPI_Wait( &rrequest, &rstat);
    }
    MPI_Buffer_detach( &b, &bufsize);
    free(b);

} // END exchangeField_movewin
