#include "cField3D.h"

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for cField3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
cField3D::cField3D() : cField()
{
    cdata_=NULL;
}

// with the dimensions as input argument
cField3D::cField3D( vector<unsigned int> dims ) : cField( dims )
{
    cdata_=NULL;
    allocateDims( dims );
}

// with the dimensions and output (dump) file name as input argument
cField3D::cField3D( vector<unsigned int> dims, string name_in ) : cField( dims, name_in )
{
    cdata_=NULL;
    allocateDims( dims );
}

// with the dimensions as input argument
cField3D::cField3D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : cField( dims, mainDim, isPrimal )
{
    cdata_=NULL;
    allocateDims( dims, mainDim, isPrimal );
}

// with the dimensions and output (dump) file name as input argument
cField3D::cField3D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name_in ) : cField( dims, mainDim, isPrimal, name_in )
{
    cdata_=NULL;
    allocateDims( dims, mainDim, isPrimal );
}

// without allocating
cField3D::cField3D( string name_in, vector<unsigned int> dims ) : cField( dims, name_in )
{
    cdata_=NULL;
    dims_=dims;
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for cField3D
// ---------------------------------------------------------------------------------------------------------------------
cField3D::~cField3D()
{

    if( cdata_!=NULL ) {
        delete [] cdata_;
        for( unsigned int i=0; i<dims_[0]; i++ ) {
            delete [] data_3D[i];
        }
        delete [] data_3D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField3D
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::allocateDims()
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if( dims_.size()!=3 ) {
        ERROR( "Alloc error must be 3 : " << dims_.size() );
    }
    if( cdata_!=NULL ) {
        delete [] cdata_;
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    cdata_ = new complex<double>[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!! (JD)}
    
    data_3D= new complex<double> **[dims_[0]];
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_3D[i]= new complex<double> *[dims_[1]];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            data_3D[i][j] = cdata_ + i*dims_[1]*dims_[2] + j*dims_[2];
            for( unsigned int k=0; k<dims_[2]; k++ ) {
                data_3D[i][j][k] = 0.0;
            }
        }
    }
    globalDims_ = dims_[0]*dims_[1]*dims_[2];
    
}

void cField3D::deallocateDataAndSetTo( Field* f )
{
    delete [] cdata_;
    cdata_ = NULL;
    delete [] data_3D;
    data_3D = NULL;
    
    cdata_ = (static_cast<cField3D *>(f))->cdata_;
    data_3D = (static_cast<cField3D *>(f))->data_3D;

}

void cField3D::allocateDims( unsigned int dims1, unsigned int dims2, unsigned int dims3 )
{
    vector<unsigned int> dims( 3 );
    dims[0]=dims1;
    dims[1]=dims2;
    dims[2]=dims3;
    allocateDims( dims );
}

// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField3D
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    if( dims_.size()!=3 ) {
        ERROR( "Alloc error must be 3 : " << dims_.size() );
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
    
    cdata_ = new complex<double>[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!! (JD)}
    
    data_3D= new complex<double> **[dims_[0]];
    for( unsigned int i=0; i<dims_[0]; i++ )  {
        data_3D[i]= new complex<double> *[dims_[1]];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            data_3D[i][j] = cdata_ + i*dims_[1]*dims_[2] + j*dims_[2];
            for( unsigned int k=0; k<dims_[2]; k++ ) {
                data_3D[i][j][k] = 0.0;
            }
        }
        
    }
    
    globalDims_ = dims_[0]*dims_[1]*dims_[2];
    
}



// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void cField3D::shift_x( unsigned int delta )
{
    memmove( &( data_3D[0][0] ), &( data_3D[delta][0] ), ( dims_[1]*dims_[0]-delta*dims_[1] )*sizeof( complex<double> ) );
    memset( &( data_3D[dims_[0]-delta][0] ), 0, delta*dims_[1]*sizeof( complex<double> ) );
    
}

double cField3D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{
    double nrj( 0. );
    
    int idxlocalstart[3];
    int idxlocalend[3];
    for( int i=0 ; i<3 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    
    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        for( int j=idxlocalstart[1] ; j<idxlocalend[1] ; j++ ) {
            for( int k=idxlocalstart[2] ; k<idxlocalend[2] ; k++ ) {
                nrj += ( data_3D[i][j][k] ).real()*( data_3D[i][j][k] ).real()+ ( data_3D[i][j][k] ).imag()*( data_3D[i][j][k] ).imag();
            }
        }
    }
    
    return nrj;
}


void cField3D::put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField3D *out3D = static_cast<cField3D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.oversize[1] ) ;
    int kout = thisPatch->Pcoordinates[2]*params.n_space[2] - ( outPatch->getCellStartingGlobalIndex(2) + params.oversize[2] ) ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        for( unsigned int j = 0 ; j < this->dims_[1] ; j++ ) {
            for( unsigned int k = 0 ; k < this->dims_[2] ; k++ ) {
                ( *out3D )( iout+i, jout+j, kout+k ) = ( *this )( i, j, k );
            }
        }
    }
    
}

void cField3D::add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField3D *out3D = static_cast<cField3D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.oversize[1] ) ;
    int kout = thisPatch->Pcoordinates[2]*params.n_space[2] - ( outPatch->getCellStartingGlobalIndex(2) + params.oversize[2] ) ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        for( unsigned int j = 0 ; j < this->dims_[1] ; j++ ) {
            for( unsigned int k = 0 ; k < this->dims_[2] ; k++ ) {
                ( *out3D )( iout+i, jout+j, kout+k ) += ( *this )( i, j, k );
            }
        }
    }
    
}

void cField3D::get( Field *inField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *thisPatch )
{
    cField3D *in3D  = static_cast<cField3D *>( inField );
    
    std::vector<unsigned int> dual =  in3D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.n_space[0] - ( inPatch->getCellStartingGlobalIndex(0) + params.oversize[0] );
    int jin = thisPatch->Pcoordinates[1]*params.n_space[1] - ( inPatch->getCellStartingGlobalIndex(1) + params.oversize[1] );
    int kin = thisPatch->Pcoordinates[2]*params.n_space[2] - ( inPatch->getCellStartingGlobalIndex(2) + params.oversize[2] );
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        for( unsigned int j = 0 ; j < this->dims_[1] ; j++ ) {
            for( unsigned int k = 0 ; k < this->dims_[2] ; k++ ) {
                ( *this )( i, j, k ) = ( *in3D )( iin+i, jin+j, kin+k );
            }
        }
    }
}
