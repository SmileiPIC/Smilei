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
void Field1D::allocateDims(std::vector<unsigned int> dims )
{
    // for Field1D only
	if (dims.size()!=1) ERROR("Alloc error must be 1 : " << dims_.size());
	
	dims_ = dims;
	data_ = new double[ dims[0] ];
	//! \todo{change to memset (JD)}
	for (unsigned int i=0;i<dims_[0];i++) data_[i]=0.0;	
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field1D
// ---------------------------------------------------------------------------------------------------------------------
void Field1D::dump(vector<unsigned int> dims)
{
  //cout << "print " << dims[0] << " elements " << endl;
	for (unsigned int i=0 ; i<dims[0] ; i++)  fdata_ << data_[i] << endl;
	fdata_ << endl;
}
