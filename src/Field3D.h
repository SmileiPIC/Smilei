#ifndef Field3D_H
#define Field3D_H

#include "Field.h"
#include <vector>



//! class Field3D used to defined a 3d vector
class Field3D : public Field
{
    
public:
    //! Constructor for Field3D: no input argument
	Field3D();
    
    //! Constructor for Field2D: with the vector dimension as input argument
	Field3D( std::vector<unsigned int> dims );
    
	//! Constructor for Field2D: with the vector dimension and filename for the dump as input argument
    Field3D( std::vector<unsigned int> dims, std::string name );
	
    //! Destructor for Field3D
    ~Field3D();

    //! Method used to allocate a Field3D
	void allocateDims(std::vector<unsigned int> dims );

    //! Method used to dump the data contained in a Field3D
	void dump(std::vector<unsigned int> dims);
    
    //! Overloading of the () operator allowing to set a new value for the (i,j,k) element of a Field3D
	inline double& operator () (unsigned int i,unsigned int j,unsigned int k)
    {
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR("Out of limits "<< i << " " << j << " " << k));
		return data_[i][j][k]; 
	};

    //! Overloading of the () operator allowing to get the value for the (i,j,k) element of a Field3D
	inline double operator () (unsigned int i,unsigned int j,unsigned int k) const { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR("Out of limits "<< i << " " << j << " " << k));
		return data_[i][j][k]; 
	};

private:
    //!\todo{Comment what are these stuffs (MG for JD)}
	double*** data_;
	double *data_3D;
};

#endif
