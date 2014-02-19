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
    data_2D=NULL;
}

// with the dimensions as input argument
Field2D::Field2D(vector<unsigned int> dims) : Field(dims)
{
    data_2D=NULL;
    allocateDims(dims);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D(vector<unsigned int> dims, string name) : Field(dims, name)
{
    data_2D=NULL;
    allocateDims(dims);
}

// with the dimensions as input argument
Field2D::Field2D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) : Field(dims, mainDim, isPrimal)
{
    data_2D=NULL;
    allocateDims(dims, mainDim, isPrimal);
}

// with the dimensions and output (dump) file name as input argument
Field2D::Field2D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name) : Field(dims, mainDim, isPrimal, name)
{
    data_2D=NULL;
    allocateDims(dims, mainDim, isPrimal);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field2D
// ---------------------------------------------------------------------------------------------------------------------
Field2D::~Field2D()
{

    if (data_2D!=NULL) {
        delete [] data_2D;
        delete [] data_;
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
    if (data_2D!=NULL) delete [] data_2D;

    isPrimal_.resize( dims.size(), 0 );

    data_2D = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}

    data_= new double*[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++) {
        data_[i] = data_2D + i*dims_[1];
        for (unsigned int j=0; j<dims_[1]; j++) data_[i][j] = 0.0;
    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
    dims_=dims;
    if (dims_.size()!=2) ERROR("Alloc error must be 2 : " << dims.size());
    if (data_2D) delete [] data_2D;

    // isPrimal define if mainDim is Primal or Dual
    // isPrimal_ = 0 if Prim  = 1 if Dual
    isPrimal_.resize( dims.size(), 0 );
    for ( unsigned int j=0 ; j<dims.size() ; j++ ) {
        if ( (j==mainDim) && (!isPrimal) )
            isPrimal_[j] = 1;
        else if ( (j!=mainDim) && (isPrimal) )
            isPrimal_[j] = 1;
    }

    for ( unsigned int j=0 ; j<dims.size() ; j++ )
        dims_[j] += isPrimal_[j];

    data_2D = new double[dims_[0]*dims_[1]];
    //! \todo{check row major order!!! (JD)}

    data_= new double*[dims_[0]];
    for (unsigned int i=0; i<dims_[0]; i++)  {
        data_[i] = data_2D + i*dims_[1];
        for (unsigned int j=0; j<dims_[1]; j++) data_[i][j] = 0.0;
    }

}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::dump(vector<unsigned int> dims)
{

}

