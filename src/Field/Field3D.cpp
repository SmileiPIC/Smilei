#include "Field3D.h"

#include <iostream>
#include <vector>
#include <cstring>

#include "Params.h"
#include "SmileiMPI.h"
#include "Patch.h"
#include "Tools.h"

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field3D::Field3D() : Field()
{
    data_=NULL;
}

// with the dimensions as input argument
Field3D::Field3D( vector<unsigned int> dims ) : Field( dims )
{
    data_=NULL;
    allocateDims( dims );
}

// with the dimensions and output (dump) file name as input argument
Field3D::Field3D( vector<unsigned int> dims, string name_in ) : Field( dims, name_in )
{
    data_=NULL;
    allocateDims( dims );
}

// with the dimensions as input argument
Field3D::Field3D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : Field( dims, mainDim, isPrimal )
{
    data_=NULL;
    allocateDims( dims, mainDim, isPrimal );
}

// with the dimensions and output (dump) file name as input argument
Field3D::Field3D( vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name_in ) : Field( dims, mainDim, isPrimal, name_in )
{
    data_=NULL;
    allocateDims( dims, mainDim, isPrimal );
}


// without allocating
Field3D::Field3D( string name_in, vector<unsigned int> dims ) : Field( dims, name_in )
{
    data_=NULL;
    dims_=dims;
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field3D
// ---------------------------------------------------------------------------------------------------------------------
Field3D::~Field3D()
{
    if( data_!=NULL ) {
        delete [] data_;
        for( unsigned int i=0; i<dims_[0]; i++ ) {
            delete [] this->data_3D[i];
        }
        delete [] this->data_3D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field3D
// ---------------------------------------------------------------------------------------------------------------------
void Field3D::allocateDims()
{
    if( dims_.size()!=3 ) {
        ERROR( "Alloc error must be 3 : " << dims_.size() );
    }
    if( data_ ) {
        delete [] data_;
    }
    
    isDual_.resize( dims_.size(), 0 );
    
    data_ = new double[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!!}
    data_3D= new double **[dims_[0]];
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_3D[i]= new double*[dims_[1]];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            data_3D[i][j] = data_ + i*dims_[1]*dims_[2] + j*dims_[2];
            for( unsigned int k=0; k<dims_[2]; k++ ) {
                this->data_3D[i][j][k] = 0.0;
            }
        }
    }//i
    
    globalDims_ = dims_[0]*dims_[1]*dims_[2];
    
}

void Field3D::deallocateDataAndSetTo( Field* f )
{
    delete [] data_;
    data_ = NULL;
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        delete [] data_3D[i];
    }
    delete [] data_3D;
    data_3D = NULL;

    data_   = f->data_;
    data_3D = (static_cast<Field3D *>(f))->data_3D;
    
}


void Field3D::allocateDims( unsigned int dims1, unsigned int dims2, unsigned int dims3 )
{
    vector<unsigned int> dims( 3 );
    dims[0]=dims1;
    dims[1]=dims2;
    dims[2]=dims3;
    allocateDims( dims );
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field3D
// ---------------------------------------------------------------------------------------------------------------------
void Field3D::allocateDims( unsigned int mainDim, bool isPrimal )
{
    if( dims_.size()!=3 ) {
        ERROR( "Alloc error must be 3 : " << dims_.size() );
    }
    if( data_ ) {
        delete [] data_;
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
    
    data_ = new double[dims_[0]*dims_[1]*dims_[2]];
    //! \todo{check row major order!!!}
    data_3D= new double **[dims_[0]*dims_[1]];
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        data_3D[i]= new double*[dims_[1]];
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            this->data_3D[i][j] = data_ + i*dims_[1]*dims_[2] + j*dims_[2];
            for( unsigned int k=0; k<dims_[2]; k++ ) {
                this->data_3D[i][j][k] = 0.0;
            }
        }
    }//i
    
    globalDims_ = dims_[0]*dims_[1]*dims_[2];
    
    //isDual_ = isPrimal;
}



// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void Field3D::shift_x( unsigned int delta )
{
    memmove( &( data_3D[0][0][0] ), &( data_3D[delta][0][0] ), ( dims_[2]*dims_[1]*dims_[0]-delta*dims_[2]*dims_[1] )*sizeof( double ) );
    memset( &( data_3D[dims_[0]-delta][0][0] ), 0, delta*dims_[1]*dims_[2]*sizeof( double ) );
    
}

double Field3D::norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] )
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
                nrj += data_3D[i][j][k]*data_3D[i][j][k];
            }
        }
    }
    
    return nrj;
}

void Field3D::extract_slice_yz( unsigned int ix, Field2D *slice )
{
    DEBUGEXEC( if( dims_[1]!=slice->dims_[0] ) ERROR( name << " : " <<  dims_[1] << " and " << slice->dims_[0] ) );
    DEBUGEXEC( if( dims_[2]!=slice->dims_[1] ) ERROR( name << " : " <<  dims_[2] << " and " << slice->dims_[1] ) );
    
    for( unsigned int j=0; j<dims_[1]; j++ ) {
        for( unsigned int k=0; k<dims_[2]; k++ ) {
            ( *slice )( j, k ) = ( *this )( ix, j, k );
        }
    }
    
}

void Field3D::extract_slice_xz( unsigned int iy, Field2D *slice )
{
    DEBUGEXEC( if( dims_[0]!=slice->dims_[0] ) ERROR( name << " : " <<  dims_[0] << " and " << slice->dims_[0] ) );
    DEBUGEXEC( if( dims_[2]!=slice->dims_[1] ) ERROR( name << " : " <<  dims_[2] << " and " << slice->dims_[1] ) );
    
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        for( unsigned int k=0; k<dims_[2]; k++ ) {
            ( *slice )( i, k ) = ( *this )( i, iy, k );
        }
    }
    
}

void Field3D::extract_slice_xy( unsigned int iz, Field2D *slice )
{
    DEBUGEXEC( if( dims_[0]!=slice->dims_[0] ) ERROR( name << " : " <<  dims_[0] << " and " << slice->dims_[0] ) );
    DEBUGEXEC( if( dims_[1]!=slice->dims_[1] ) ERROR( name << " : " <<  dims_[1] << " and " << slice->dims_[1] ) );
    
    for( unsigned int i=0; i<dims_[0]; i++ ) {
        for( unsigned int j=0; j<dims_[1]; j++ ) {
            ( *slice )( i, j ) = ( *this )( i, j, iz );
        }
    }
    
}


void Field3D::put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    Field3D *out3D = static_cast<Field3D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    int kout = thisPatch->Pcoordinates[2]*params.n_space[2] - ( outPatch->getCellStartingGlobalIndex(2) + params.region_oversize[2] ) ;    
    
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            for( unsigned int k = 0 ; k < params.n_space[2]+1+dual[2]+2*params.oversize[2] ; k++ ) {
                ( *out3D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1], kout+k+params.region_oversize[2]-params.oversize[2] ) = ( *this )( i, j, k );
            }
        }
    }
    
}


void Field3D::add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch )
{
    Field3D *out3D = static_cast<Field3D *>( outField );
    
    std::vector<unsigned int> dual =  this->isDual_;
    
    int iout = thisPatch->Pcoordinates[0]*params.n_space[0] - ( outPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] ) ;
    int jout = thisPatch->Pcoordinates[1]*params.n_space[1] - ( outPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] ) ;
    int kout = thisPatch->Pcoordinates[2]*params.n_space[2] - ( outPatch->getCellStartingGlobalIndex(2) + params.region_oversize[2] ) ;
    
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            for( unsigned int k = 0 ; k < params.n_space[2]+1+dual[2]+2*params.oversize[2] ; k++ ) {
                ( *out3D )( iout+i+params.region_oversize[0]-params.oversize[0], jout+j+params.region_oversize[1]-params.oversize[1], kout+k+params.region_oversize[2]-params.oversize[2] ) += ( *this )( i, j, k );
            }
        }
    }
    
}

void Field3D::get( Field *inField, Params &params, SmileiMPI *smpi, Patch *inPatch, Patch *thisPatch )
{
    Field3D *in3D  = static_cast<Field3D *>( inField );
    
    std::vector<unsigned int> dual =  in3D->isDual_;
    
    int iin = thisPatch->Pcoordinates[0]*params.n_space[0] - ( inPatch->getCellStartingGlobalIndex(0) + params.region_oversize[0] );
    int jin = thisPatch->Pcoordinates[1]*params.n_space[1] - ( inPatch->getCellStartingGlobalIndex(1) + params.region_oversize[1] );
    int kin = thisPatch->Pcoordinates[2]*params.n_space[2] - ( inPatch->getCellStartingGlobalIndex(2) + params.region_oversize[2] );
    
    for( unsigned int i = 0 ; i < params.n_space[0]+1+dual[0]+2*params.oversize[0] ; i++ ) {
        for( unsigned int j = 0 ; j < params.n_space[1]+1+dual[1]+2*params.oversize[1] ; j++ ) {
            for( unsigned int k = 0 ; k < params.n_space[2]+1+dual[2]+2*params.oversize[2] ; k++ ) {
                ( *this )( i, j, k ) = ( *in3D )( iin+i+params.region_oversize[0]-params.oversize[0], jin+j+params.region_oversize[1]-params.oversize[1], kin+k+params.region_oversize[2]-params.oversize[2] );
            }
        }
    }
    
}
