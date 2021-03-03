
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
                    ntype_[ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
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
                ntype_[ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
            }
        }
    }
    
}


Patch3D::~Patch3D()
{
    cleanType();
}

// ---------------------------------------------------------------------------------------------------------------------
// Create MPI_Datatypes used in initSumField and initExchange
// ---------------------------------------------------------------------------------------------------------------------
void Patch3D::createType2( Params &params )
{
    //int nx0 = params.n_space_region[0] + 1 + 2*oversize[0];
    int ny0 = params.n_space_region[1] + 1 + 2*oversize[1];
    int nz0 = params.n_space_region[2] + 1 + 2*oversize[2];
    
    int ny, nz;
    //int nx_sum, ny_sum, nz_sum;
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        //nx = nx0 + ix_isPrim;
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            ny = ny0 + iy_isPrim;
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                nz = nz0 + iz_isPrim;
                
                // Still used ???
                ntype_[ix_isPrim][iy_isPrim][iz_isPrim] = MPI_DATATYPE_NULL;
                if (!params.is_pxr)
                    MPI_Type_contiguous(nz*ny*params.n_space[0], MPI_DOUBLE, &(ntype_[ix_isPrim][iy_isPrim][iz_isPrim]));   //clrw lines
                else
                    MPI_Type_contiguous(nz*ny*(params.n_space[0]), MPI_DOUBLE, &(ntype_[ix_isPrim][iy_isPrim][iz_isPrim]));   //clrw lines
                MPI_Type_commit( &(ntype_[ix_isPrim][iy_isPrim][iz_isPrim]) );
            }
        }
    }
    
}

void Patch3D::cleanType()
{
    if( ntype_[0][0][0] == MPI_DATATYPE_NULL ) {
        return;
    }
    
    for( int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++ ) {
        for( int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++ ) {
            for( int iz_isPrim=0 ; iz_isPrim<2 ; iz_isPrim++ ) {
                MPI_Type_free( &(ntype_[ix_isPrim][iy_isPrim][iz_isPrim]) );
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

    MPI_Datatype ntype = ntype_[isDual[0]][isDual[1]][isDual[2]]; //ntype_[3] is clrw plans.
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
