#include "Field1D.h"

#include <iostream>
#include <vector>
#include <cstring>

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field1D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field1D::Field1D() : Field()
{
}

// with the dimensions as input argument
Field1D::Field1D(vector<unsigned int> dims) : Field(dims)
{
    allocateDims(dims);
}

// with the dimensions and output (dump) file name as input argument
Field1D::Field1D(vector<unsigned int> dims, string name) : Field(dims, name)
{
    allocateDims(dims);
}

// with the dimensions as input argument
Field1D::Field1D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) : Field(dims, mainDim, isPrimal)
{
    allocateDims(dims, mainDim, isPrimal);
}

// with the dimensions and output (dump) file name as input argument
Field1D::Field1D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name) : Field(dims, mainDim, isPrimal, name)
{
    allocateDims(dims, mainDim, isPrimal);
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field1D
// ---------------------------------------------------------------------------------------------------------------------
Field1D::~Field1D()
{
    if (data_!=NULL) {
        delete [] data_;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::allocateDims(std::vector<unsigned int> dims)
{
    // for Field1D only
    dims_ = dims;
    if (dims.size()!=1) ERROR("Alloc error must be 1 : " << dims.size());

    isDual_.resize( dims.size(), 0 );

    data_ = new double[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for (unsigned int i=0; i<dims_[0]; i++) data_[i]=0.0;

    globalDims_ = dims_[0];

}

void Field1D::allocateDims(unsigned int dims1)
{
	vector<unsigned int> dims(1);
	dims[0]=dims1;
	allocateDims(dims);
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal)
{
    // for Field1D only
    dims_ = dims;
    if (dims.size()!=1) ERROR("Alloc error must be 1 : " << dims.size());

    // isPrimal define if mainDim is Primal or Dual
    isDual_.resize( dims.size(), 0 );
    for ( unsigned int j=0 ; j<dims.size() ; j++ ) {
        if ( (j==mainDim) && (!isPrimal) )
            isDual_[j] = 1;
        else if ( (j!=mainDim) && (isPrimal) )
            isDual_[j] = 1;
    }

    for ( unsigned int j=0 ; j<dims.size() ; j++ )
        dims_[j] += isDual_[j];

    data_ = new double[ dims_[0] ];
    //! \todo{change to memset (JD)}
    for (unsigned int i=0; i<dims_[0]; i++) data_[i]=0.0;

    globalDims_ = dims_[0];

}


// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::dump(vector<unsigned int> dims)
{
}


// ---------------------------------------------------------------------------------------------------------------------
// Method to shift field in space
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::shift_x(unsigned int delta)
{
    memmove( &(data_[0]), &(data_[delta]), (dims_[0]-delta)*sizeof(double) );
    //memset ( &(data_[dims_[0]-delta]), 0, delta*sizeof(double));
    for (int i=dims_[0]-delta;i<(int)dims_[0];i++) data_[i] = 0.;

}

// ---------------------------------------------------------------------------------------------------------------------
// Method used to compute the field energy (assuming electromagnetic type)
// ---------------------------------------------------------------------------------------------------------------------
double Field1D::computeNRJ(unsigned int shift, unsigned int istart[3][2], unsigned int bufsize[3][2]) {
    
    double nrj(0.);
    
    int idxlocalstart[1];
    int idxlocalend[1];
    for ( int i=0 ; i<1 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    idxlocalend[0] = istart[0][isDual_[0]]+shift;
    
    if (!isDual_[0]) {
        // the total energy is just the square of all fields
        for ( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
            nrj += data_[i]*data_[i];
        }
    } else {
        // the total energy is more complex to compute
        nrj += 0.5*data_[idxlocalend[0]-1];
        for ( int i=idxlocalstart[0] ; i<idxlocalend[0]-1 ; i++ ) {
            nrj += data_[i]*data_[i] + data_[i+1]*data_[i];
        }
    }
    
    nrj = 1.0;
    
    
    /*
    int idxlocalstart[1];
    int idxlocalend[1];
    for ( int i=0 ; i<1 ; i++ ) {
        idxlocalstart[i] = istart[i][isDual_[i]];
        idxlocalend[i]   = istart[i][isDual_[i]]+bufsize[i][isDual_[i]];
    }
    idxlocalend[0] = istart[0][isDual_[0]]+shift;
    
    for ( int i=idxlocalstart[0] ; i<idxlocalend[0] ; i++ ) {
        nrj += data_[i]*data_[i];
    }
     */
    
    return nrj;
}
