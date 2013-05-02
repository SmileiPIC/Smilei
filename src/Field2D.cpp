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
Field2D::Field2D(vector<unsigned int> dims, string name) : Field(dims)
{
	allocateDims(dims);
	fdata_.open(name.c_str(), ios::out | ios::binary);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field2D
// ---------------------------------------------------------------------------------------------------------------------
Field2D::~Field2D()
{
	delete [] data_2D;
	
	DEBUG(10,"Field 2D deleted");
	fdata_.close();
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::allocateDims(std::vector<unsigned int> dims )
{
    //! \todo{Comment on what you are doing here (MG for JD)}
	dims_=dims;
	if (dims_.size()!=2) ERROR("Alloc error must be 2 : " << dims_.size());
	if (data_2D) delete [] data_2D;

	data_2D = new double[dims_[0]*dims_[1]];
	//! \todo{check row major order!!! (JD)}
    
	data_= new double*[dims_[0]];
	for (unsigned int i=0; i<dims_[0]; i++)
    {
		data_[i] = data_2D + i*dims_[1];
	}
	
	DEBUG(10,"Fields 2D created: " << data_[0] << "x" << data_[1]);
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field2D
// ---------------------------------------------------------------------------------------------------------------------
void Field2D::dump(vector<unsigned int> dims)
{
	for (unsigned int i=0 ; i<dims[0] ; i++)    fdata_ << data_[i] << " ";
	fdata_ << endl;
    
}

