#include "Field3D.h"

#include <iostream>
#include <vector>

using namespace std;



// ---------------------------------------------------------------------------------------------------------------------
// Creators for Field3D
// ---------------------------------------------------------------------------------------------------------------------

// with no input argument
Field3D::Field3D() : Field()
{
	data_3D=NULL;
}

// with the dimensions as input argument
Field3D::Field3D(vector<unsigned int> dims) : Field(dims)
{
	data_3D=NULL;
	allocateDims(dims);
}

// with the dimensions and output (dump) file name as input argument
Field3D::Field3D(vector<unsigned int> dims, string name) : Field(dims)
{
	allocateDims(dims);
	fdata_.open(name.c_str(), ios::out | ios::binary);
}



// ---------------------------------------------------------------------------------------------------------------------
// Destructor for Field3D
// ---------------------------------------------------------------------------------------------------------------------
Field3D::~Field3D()
{
	delete [] data_3D;
	for (unsigned int i=0; i<dims_[0]; i++) delete data_[i];
	delete data_;
	
	fdata_.close();
	//DEBUG(10,"Field 3D deleted");
}



// ---------------------------------------------------------------------------------------------------------------------
// Method used for allocating the dimension of a Field3D
// ---------------------------------------------------------------------------------------------------------------------
void Field3D::allocateDims(std::vector<unsigned int> dims ) {
	dims_=dims;
	if (dims_.size()!=3) ERROR("Alloc error must be 3 : " << dims_.size());

	if (data_3D) delete [] data_3D;

	data_3D = new double[dims_[0]*dims_[1]*dims_[2]];
	//! \todo{check row major order!!!}
	data_= new double**[dims_[0]*dims_[1]];
	for (unsigned int i=0; i<dims_[0]; i++)
    {
		data_[i]= new double*[dims_[1]];
		for (unsigned int j=0; j<dims_[1]; j++)
        {
			data_[i][j] = data_3D + i*dims_[1]*dims_[2] + j*dims_[2];
		}
	}//i
	
	//DEBUG(10,"Fields 3D created: " << dims_[0] << "x" << dims_[1] << "x" << dims_[2]);

}



// ---------------------------------------------------------------------------------------------------------------------
// Method used to create a dump for a Field3D
// ---------------------------------------------------------------------------------------------------------------------
void Field3D::dump(vector<unsigned int> dims)
{
	for (unsigned int i=0 ; i<dims[0] ; i++)    fdata_ << data_[i] << " ";
	fdata_ << endl;
    
}
