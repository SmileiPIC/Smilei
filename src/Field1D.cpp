#include "Field1D.h"

#include <iostream>
#include <vector>
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
Field1D::Field1D(vector<unsigned int> dims, string name) : Field(dims)
{
	allocateDims(dims);
	fdata_.open(name.c_str(), ios::out);
}

// with the dimensions as input argument
Field1D::Field1D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) : Field(dims, mainDim, isPrimal)
{
	allocateDims(dims, mainDim, isPrimal);
}

// with the dimensions and output (dump) file name as input argument
Field1D::Field1D(vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, string name) : Field(dims, mainDim, isPrimal)
{
	allocateDims(dims, mainDim, isPrimal);
	fdata_.open(name.c_str(), ios::out);
}


// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field1D
// ---------------------------------------------------------------------------------------------------------------------
Field1D::~Field1D()
{
	delete [] data_;
	fdata_.close();
}


// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::allocateDims(std::vector<unsigned int> dims)
{
	// for Field1D only
	dims_ = dims;
	if (dims.size()!=1) ERROR("Alloc error must be 1 : " << dims.size());

	isPrimal_.resize( dims.size(), 0 );

	data_ = new double[ dims_[0] ];
	//! \todo{change to memset (JD)}
	for (unsigned int i=0;i<dims_[0];i++) data_[i]=0.0;

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

	data_ = new double[ dims_[0] ];
	//! \todo{change to memset (JD)}
	for (unsigned int i=0;i<dims_[0];i++) data_[i]=0.0;

}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::dump(vector<unsigned int> dims)
{
  fdata_.precision( 20 ); // floatfield set to fixed
  //cout << "print " << dims[0] << " elements " << endl;
	for (unsigned int i=0 ; i<dims[0] ; i++)  fdata_ << data_[i] << endl;
	fdata_ << endl;
}
