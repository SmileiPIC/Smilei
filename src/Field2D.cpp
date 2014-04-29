#include "Field2D.h"

#include <iostream>
#include <vector>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field2D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field2D::Field2D() : Field()
{
    data_=NULL;
}

// with the dimensions as input argument
Field2D::Field2D(vector<unsigned int> dims) : Field(dims)
{
    data_=NULL;
    allocateDims(dims);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D(vector<unsigned int> dims, string name) : Field(dims, name)
{
    data_=NULL;
    allocateDims(dims);
}

// with the dimensions as input argument
Field2D::Field2D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) : Field(dims, mainDim, isPrimal)
{
    data_=NULL;
    allocateDims(dims, mainDim, isPrimal);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name) : Field(dims, mainDim, isPrimal, name)
{
    data_=NULL;
    allocateDims(dims, mainDim, isPrimal);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field2D
// ---------------------------------------------------------------------------------------------------------------------
Field2D::~Field2D()
{

    if (data_!=NULL) {
        delete [] data_;
        delete [] data_2D;
    }
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims(std::vector<unsigned int> dims )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    dims_=dims;
    if (dims_.size()!=2) ERROR("Alloc error must be 2 : " << dims.size());
    if (data_!=NULL) delete [] data_;
	
    isDual_.resize( dims.size(), 0 );
	
    data_ = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}
	
    data_2D= new double*[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++) {
        data_2D[i] = data_ + i*dims_[1];
        for (unsigned int j=0; j<dims_[1]; j++) data_2D[i][j] = 0.0;
    }
	
    globalDims_ = dims_[0]*dims_[1];
	
}

void Field2D::allocateDims(unsigned int dims1, unsigned int dims2)
{
	vector<unsigned int> dims(2);
	dims[0]=dims1;
	dims[1]=dims2;
	allocateDims(dims);
}

// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    dims_=dims;
    if (dims_.size()!=2) ERROR("Alloc error must be 2 : " << dims.size());
    if (data_) delete [] data_;

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

    data_ = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}

    data_2D= new double*[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++)  {
        data_2D[i] = data_ + i*dims_[1];
        for (unsigned int j=0; j<dims_[1]; j++) data_2D[i][j] = 0.0;
    }

    globalDims_ = dims_[0]*dims_[1];

}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::dump(vector<unsigned int> dims)
{

}

