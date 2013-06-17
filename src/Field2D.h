#ifndef Field2D_H
#define Field2D_H

#include "Field.h"
#include <vector>



//! class Field2D used to defined a 2d vector
class Field2D : public Field
{
    
 public:
	//! Constructor for Field2D: no input argument
	Field2D();
    
	//! Constructor for Field2D: with the vector dimension as input argument
	Field2D( std::vector<unsigned int> dims );
	Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
	//! Constructor for Field2D: with the vector dimension and filename for the dump as input argument
	Field2D( std::vector<unsigned int> dims, std::string name );
	Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );
    
	//! Destructor for Field2D
	~Field2D();

	//! Method used to allocate a Field2D
	void allocateDims(std::vector<unsigned int> dims );
	void allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
	//! Method used to dump the data contained in a Field2D
	void dump(std::vector<unsigned int> dims);

	//! Overloading of the () operator allowing to set a new value for the (i,j) element of a Field2D
	inline double& operator () (unsigned int i,unsigned int j) {
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR("Out of limits "<< i << " " << j));
		return data_[i][j]; 
	};

	//! Overloading of the () operator allowing to get the value of the (i,j) element of a Field2D
	inline double operator () (unsigned int i,unsigned int j) const { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR("Out of limits "<< i << " " << j));
		return data_[i][j]; 
	};
	double** data_;
	//! \todo{for debbugging, to remove (JD)}
	inline void setData_(int i, double val) {data_2D[i] = val;}
	
 private:
	//!\todo{Comment what are these stuffs (MG for JD)}
	double *data_2D;
};

#endif

