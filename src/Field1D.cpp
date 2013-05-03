#include "Field1D.h"

#include <iostream>
#include <vector>
using namespace std;

Field1D::Field1D() : Field()
{
}


Field1D::Field1D(vector<unsigned int> dims) : Field(dims)
{
	allocateDims(dims);
}


Field1D::Field1D(vector<unsigned int> dims, string name) : Field(dims)
{
	allocateDims(dims);
	
	fdata_.open(name.c_str(), ios::out);
	
}

void Field1D::allocateDims(std::vector<unsigned int> dims ) {
	if (dims.size()!=1) ERROR("Alloc error must be 1 : " << dims_.size());
	
	dims_ = dims;
	data_ = new double[ dims[0] ];
	//! \todo{change to memset}
	for (unsigned int i=0;i<dims_[0];i++) data_[i]=0.0;	
}


Field1D::~Field1D()
{
	delete [] data_;
	fdata_.close();
}


void Field1D::dump(vector<unsigned int> dims)
{
	for (unsigned int i=0 ; i<dims[0] ; i++)
		fdata_ << data_[i] << endl;
	fdata_ << endl;
	fdata_ << endl;
}
