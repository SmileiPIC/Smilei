
#include "Field2D.h"

#include <iostream>
#include <vector>
using namespace std;

Field2D::Field2D() : Field()
{
	data_2D=NULL;
}

Field2D::Field2D(vector<unsigned int> dims) : Field(dims)
{
	data_2D=NULL;
	allocateDims(dims);
}

void Field2D::allocateDims(std::vector<unsigned int> dims ) {
	dims_=dims;
	if (dims_.size()!=2) ERROR("Alloc error must be 2 : " << dims_.size());
	if (data_2D) delete [] data_2D;

	data_2D = new double[dims_[0]*dims_[1]];
	//! \todo{check row major order!!!}
	data_= new double*[dims_[0]];
	for (unsigned int i=0; i<dims_[0]; i++) {
		data_[i] = data_2D + i*dims_[1];
	}
	
	DEBUG(10,"Fields 2D created: " << data_[0] << "x" << data_[1]);
}


Field2D::Field2D(vector<unsigned int> dims, string name) : Field(dims)
{
	allocateDims(dims);
	
	fdata_.open(name.c_str(), ios::out | ios::binary);
}


Field2D::~Field2D()
{
	delete [] data_2D;
	
	DEBUG(10,"Field 2D deleted");
	fdata_.close();
}


void Field2D::dump(vector<unsigned int> dims)
{
	for (unsigned int i=0 ; i<dims[0] ; i++)
		fdata_ << data_[i] << " ";
	fdata_ << endl;
    
}

