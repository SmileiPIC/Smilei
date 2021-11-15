#include "cField2D.h"

#include <iostream>
#include <vector>
#include <cstring>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for cField2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
cField2D::cField2D() : cField()
{
    cleaned_ = false;
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions as input argument
cField2D::cField2D( vector<unsigned int> dims ) : cField( dims )
{
    cleaned_ = false;
    allocateDims( dims );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions and output (dump) file name as input argument
cField2D::cField2D( vector<unsigned int> dims, string name_in ) : cField( dims, name_in )
{
    cleaned_ = false;
    allocateDims( dims );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions as input argument
cField2D::cField2D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : cField( dims, mainDim, isPrimal )
{
    cleaned_ = false;
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// with the dimensions and output (dump) file name as input argument
cField2D::cField2D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name_in ) : cField( dims, mainDim, isPrimal, name_in )
{
    cleaned_ = false;
    allocateDims( dims, mainDim, isPrimal );
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}

// without allocating
cField2D::cField2D( string name_in, vector<unsigned int> dims ) : cField( dims, name_in )
{
    cleaned_ = false;
    dims_ = dims;
    globalDims_ = dims_[0]*dims_[1];
    sendFields_.resize(4,NULL);
    recvFields_.resize(4,NULL);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for cField2D
// ---------------------------------------------------------------------------------------------------------------------
cField2D::~cField2D()
{
    if (cleaned_ )
        return;
    for (int iside=0 ; iside<(int)(sendFields_.size()) ; iside++ ) {
        if ( sendFields_[iside] != NULL ) {
            delete sendFields_[iside];
            sendFields_[iside] = NULL;
            delete recvFields_[iside];
            recvFields_[iside] = NULL;
        }
    }
    if( cdata_!=NULL ) {
        delete [] cdata_;
        delete [] data_2D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField2D
// ---------------------------------------------------------------------------------------------------------------------
void cField2D::allocateDims()
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if( dims_.size()!=2 ) {
        ERROR( "Alloc error must be 2 : " << dims_.size() );
    }
    if( cdata_!=NULL ) {
        delete [] cdata_;
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    cdata_ = new complex<double>[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}
    
    data_2D= new complex<double> *[dims_[0]];
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_2D[i] = cdata_ + i*dims_[1];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            data_2D[i][j] = 0.0;
        }
    }
    
    globalDims_ = dims_[0]*dims_[1];
    
}

void cField2D::deallocateDataAndSetTo( Field* f )
{
    delete [] cdata_;
    cdata_ = NULL;
    delete [] data_2D;
    data_2D = NULL;
    cleaned_ = true;

    cdata_ = (static_cast<cField2D *>(f))->cdata_;
    data_2D = (static_cast<cField2D *>(f))->data_2D;
    
}

void cField2D::allocateDims( unsigned int dims1, unsigned int dims2 )
{
    vector<unsigned int> dims( 2 );
    dims[0]=dims1;
    dims[1]=dims2;
    allocateDims( dims );
}

// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField2D
// ---------------------------------------------------------------------------------------------------------------------
void cField2D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if( dims_.size()!=2 ) {
        ERROR( "Alloc error must be 2 : " << dims_.size() );
    }
    if( cdata_ ) {
        delete [] cdata_;
    }
    
    // isPrimal define if mainDim is Primal or Dual
    isDual_.resize( dims_.size(), 0 );
    for( unsigned int j=0 ; j<dims_.size() ; j++ ) {
        if( ( j==mainDim ) && ( !isPrimal ) ) {
            isDual_[j] = 1;
        } else if( ( j!=mainDim ) && ( isPrimal ) ) {
            isDual_[j] = 1;
        }
    }
    
    for( unsigned int j=0 ; j<dims_.size() ; j++ ) {
        dims_[j] += isDual_[j];
    }
    
    cdata_ = new complex<double>[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}
    
    data_2D= new complex<double> *[dims_[0]];
    for( unsigned int i=0; i<dims_[0]; i++ )  {
        data_2D[i] = cdata_ + i*dims_[1];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            data_2D[i][j] = 0.0;
        }
    }
    
    globalDims_ = dims_[0]*dims_[1];
    
}



// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void cField2D::shift_x( unsigned int delta )
{
    memmove( &( data_2D[0][0] ), &( data_2D[delta][0] ), ( dims_[1]*dims_[0]-delta*dims_[1] )*sizeof( complex<double> ) );
    memset( &( data_2D[dims_[0]-delta][0] ), 0, delta*dims_[1]*sizeof( complex<double> ) );
    
}

double cField2D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{
    double nrj( 0. );
    
    int idxlocalstart[2];
    int idxlocalend[2];
    for( int i=0 ; i<2 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    
    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        for( int j=idxlocalstart[1] ; j<idxlocalend[1] ; j++ ) {
            nrj += ( data_2D[i][j] ).real()*( data_2D[i][j] ).real()+ ( data_2D[i][j] ).imag()*( data_2D[i][j] ).imag();
        }
    }
    
    return nrj;
}

double cField2D::norm2_cylindrical( unsigned int istart[3][2], unsigned int bufsize[3][2], int j_ref )
{
    double nrj( 0. );
    
    int idxlocalstart[2];
    int idxlocalend[2];
    for( int i=0 ; i<2 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    
    // Special case: cells on axis
    if( j_ref + idxlocalstart[1] < 0 ) {
        ERROR("IMPOSSIBLE");
    } else if( j_ref + idxlocalstart[1] == 0 ) {
        if( ! isDual_[1] ) {
            double sum = 0.;
            for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
                unsigned int j = idxlocalstart[1];
                sum += ( data_2D[i][j] ).real()*( data_2D[i][j] ).real()+ ( data_2D[i][j] ).imag()*( data_2D[i][j] ).imag();
            }
            nrj *= sum / 8.; // volume factor for on-axis cells is 1./8.
        }
        idxlocalstart[1]++;
    }
    // Remaining cells
    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        for( int j=idxlocalstart[1] ; j<idxlocalend[1] ; j++ ) {
            double volume_factor = (double)(j_ref + j) - 0.5 * isDual_[1];
            nrj += volume_factor *( ( data_2D[i][j] ).real()*( data_2D[i][j] ).real()+ ( data_2D[i][j] ).imag()*( data_2D[i][j] ).imag() );
        }
    }
    
    return nrj;
}

void cField2D::put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField2D *out2D = static_cast<cField2D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    
    //for ( unsigned int i = params.oversize[0] ; i < this->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < this->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *out2D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1] ) = ( *this )( i, j );
        }
    }
    
}

void cField2D::add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField2D *out2D = static_cast<cField2D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    
    //for ( unsigned int i = params.oversize[0] ; i < this->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < this->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *out2D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1] ) += ( *this )( i, j );
        }
    }
    
}


void cField2D::get( Field *inField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *thisPatch )
{
    cField2D *in2D  = static_cast<cField2D *>( inField );
    
    std::vector<unsigned int> dual =  in2D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.n_space[0] - ( inPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] );
    int jin = thisPatch->Pcoordinates[1]*params.n_space[1] - ( inPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] );
    
    //for ( unsigned int i = params.oversize[0] ; i < out2D->dims_[0]-params.oversize[0] ; i++ ) {
    //    for ( unsigned int j = params.oversize[1] ; j < out2D->dims_[1]-params.oversize[1] ; j++ ) {
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            ( *this )( i, j ) = ( *in2D )( iin+i+params.region_oversize[0]-params.oversize[0], jin+j+params.region_oversize[1]-params.oversize[1] );
            //( *out2D )( i, j ) = in2D->hindex;
        }
    }
    
}

void cField2D::create_sub_fields  ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> n_space = dims_;
    n_space[iDim] = ghost_size;
    if ( sendFields_[iDim*2+iNeighbor] == NULL ) {
        sendFields_[iDim*2+iNeighbor] = new cField2D(n_space);
        recvFields_[iDim*2+iNeighbor] = new cField2D(n_space);
    }
    else if ( ghost_size != int(sendFields_[iDim*2+iNeighbor]->dims_[iDim]) ) {
        delete sendFields_[iDim*2+iNeighbor];
        sendFields_[iDim*2+iNeighbor] = new cField2D(n_space);
        delete recvFields_[iDim*2+iNeighbor];
        recvFields_[iDim*2+iNeighbor] = new cField2D(n_space);
    }
}

void cField2D::extract_fields_exch( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> n_space = dims_;
    n_space[iDim] = ghost_size;

    vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * ( ghost_size + 1 + isDual_[iDim] );
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    int NX = n_space[0];
    int NY = n_space[1];

    int dimY = dims_[1];

    complex<double>* sub = static_cast<cField*>(sendFields_[iDim*2+iNeighbor])->cdata_;
    complex<double>* field = cdata_;
    for( unsigned int i=0; i< (unsigned int)(NX); i++ ) {
        for( unsigned int j=0; j<(unsigned int)(NY); j++ ) {
            sub[i*NY+j] = field[ (ix+i)*dimY+(iy+j) ];
        }
    }
}

void cField2D::inject_fields_exch ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> n_space = dims_;
    n_space[iDim] = ghost_size;

    vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - 1- ( ghost_size-1 ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    int NX = n_space[0];
    int NY = n_space[1];

    int dimY = dims_[1];

    complex<double>* sub = static_cast<cField*>(recvFields_[iDim*2+(iNeighbor+1)%2])->cdata_;
    complex<double>* field = cdata_;
    for( unsigned int i=0; i<(unsigned int)NX; i++ ) {
        for( unsigned int j=0; j<(unsigned int)NY; j++ ) {
            field[ (ix+i)*dimY+(iy+j) ] = sub[i*NY+j];
        }
    }
}

void cField2D::extract_fields_sum ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> n_space = dims_;
    n_space[iDim] = 2*ghost_size+1+isDual_[iDim];

    vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = iNeighbor * ( dims_[iDim]- ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-iNeighbor ) * 0;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    int NX = n_space[0];
    int NY = n_space[1];

    int dimY = dims_[1];

    complex<double>* sub = static_cast<cField*>(sendFields_[iDim*2+iNeighbor])->cdata_;
    complex<double>* field = cdata_;
    for( unsigned int i=0; i<(unsigned int)NX; i++ ) {
        for( unsigned int j=0; j<(unsigned int)NY; j++ ) {
            sub[i*NY+j] = field[ (ix+i)*dimY+(iy+j) ];
        }
    }
}

void cField2D::inject_fields_sum  ( int iDim, int iNeighbor, int ghost_size )
{
    std::vector<unsigned int> n_space = dims_;
    n_space[iDim] = 2*ghost_size+1+isDual_[iDim];

    vector<int> idx( 2, 0 );
    idx[iDim] = 1;
    int istart = ( ( iNeighbor+1 )%2 ) * ( dims_[iDim] - ( 2*ghost_size+1+isDual_[iDim] ) ) + ( 1-( iNeighbor+1 )%2 ) * ( 0 )  ;
    int ix = idx[0]*istart;
    int iy = idx[1]*istart;

    int NX = n_space[0];
    int NY = n_space[1];

    int dimY = dims_[1];

    complex<double>* sub = static_cast<cField*>(recvFields_[iDim*2+(iNeighbor+1)%2])->cdata_;
    complex<double>* field = cdata_;
    for( unsigned int i=0; i<(unsigned int)NX; i++ ) {
        for( unsigned int j=0; j<(unsigned int)NY; j++ ) {
            field[ (ix+i)*dimY+(iy+j) ] += sub[i*NY+j];
        }
    }
}
