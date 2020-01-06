#include "cField1D.h"

#include <iostream>
#include <vector>
#include <cstring>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creators for cField1D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
cField1D::cField1D() : cField()
{
}

// with the dimensions as input argument
cField1D::cField1D( vector<unsigned int> dims ) : cField( dims )
{
    allocateDims( dims );
}

// with the dimensions and output (dump) file name as input argument
cField1D::cField1D( vector<unsigned int> dims, string name ) : cField( dims, name )
{
    allocateDims( dims );
}

// with the dimensions as input argument
cField1D::cField1D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : cField( dims, mainDim, isPrimal )
{
    allocateDims( dims, mainDim, isPrimal );
}

// with the dimensions and output (dump) file name as input argument
cField1D::cField1D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name ) : cField( dims, mainDim, isPrimal, name )
{
    allocateDims( dims, mainDim, isPrimal );
}

// without allocating
cField1D::cField1D( string name, vector<unsigned int> dims ) : cField( dims, name )
{
    cdata_=NULL;
    dims_ = dims;
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for cField1D
// ---------------------------------------------------------------------------------------------------------------------
cField1D::~cField1D()
{
    if( cdata_!=NULL ) {
        delete [] cdata_;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField1D
// ---------------------------------------------------------------------------------------------------------------------
void cField1D::allocateDims()
{
    // for cField1D only
    if( dims_.size()!=1 ) {
        ERROR( "Alloc error must be 1 : " << dims_.size() );
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    cdata_ = new complex<double>[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        cdata_[i]=0.0;
    }
    
    globalDims_ = dims_[0];
    
}

void cField1D::deallocateDataAndSetTo( Field* f )
{
    delete [] cdata_;
    cdata_=NULL;

    cdata_ = (static_cast<cField *>(f))->cdata_;

}


void cField1D::allocateDims( unsigned int dims1 )
{
    vector<unsigned int> dims( 1 );
    dims[0]=dims1;
    allocateDims( dims );
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a cField1D
// ---------------------------------------------------------------------------------------------------------------------
void cField1D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    // for cField1D only
    if( dims_.size()!=1 ) {
        ERROR( "Alloc error must be 1 : " << dims_.size() );
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
    
    cdata_ = new complex<double>[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        cdata_[i]=0.0;
    }
    
    globalDims_ = dims_[0];
    
}


// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void cField1D::shift_x( unsigned int delta )
{
    memmove( &( cdata_[0] ), &( cdata_[delta] ), ( dims_[0]-delta )*sizeof( complex<double> ) );
    //memset ( &(data_[dims_[0]-delta]), 0, delta*sizeof(double));
    for( int i=dims_[0]-delta; i<( int )dims_[0]; i++ ) {
        cdata_[i] = 0.;
    }
    
}

double cField1D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
{

    double nrj( 0. );
    
    int idxlocalstart[1];
    int idxlocalend[1];
    for( int i=0 ; i<1 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    for( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        nrj += ( ( cdata_[i] ).real() )*( ( cdata_[i] ).real() )+( ( cdata_[i] ).imag() )*( ( cdata_[i] ).imag() );
    }
    
    return nrj;
}


void cField1D::put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField1D *out1D = static_cast<cField1D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - outPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *out1D )( iout+i ) = ( *this )( i );
    }
    
}

void cField1D::add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    cField1D *out1D = static_cast<cField1D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - outPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *out1D )( iout+i ) += ( *this )( i );
    }
    
}


void cField1D::get( Field *inField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *thisPatch )
{
    cField1D *in1D  = static_cast<cField1D *>( inField );
    
    std::vector<unsigned int> dual =  in1D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.n_space[0] - inPatch->Pcoordinates[0]*params.n_space[0]*params.global_factor[0] ;
    
    for( unsigned int i = 0 ; i < this->dims_[0] ; i++ ) {
        ( *this )( i ) = ( *in1D )( iin+i );
    }
    
}
