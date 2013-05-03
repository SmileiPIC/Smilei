
#ifndef Field2D_H
#define Field2D_H

#include "Field.h"
#include <vector>

//class Field2Dproxy {
//public:
//	friend class Field2D;
//    Field2Dproxy( double &val) : mPtr( & val ) {}
//    void operator = ( double d ) {
//    	*mPtr = d;
//    }
//    double * mPtr;
//};

class Field2D : public Field {
public:
	Field2D();
	Field2D( std::vector<unsigned int> dims );
	Field2D( std::vector<unsigned int> dims, std::string name );
	~Field2D();

	void allocateDims(std::vector<unsigned int> dims );

 	inline double& operator () (unsigned int i,unsigned int j) { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR("Out of limits "<< i << " " << j));
		return data_[i][j]; 
	};

 	inline double operator () (unsigned int i,unsigned int j) const { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR("Out of limits "<< i << " " << j));
		return data_[i][j]; 
	};
	
	void dump(std::vector<unsigned int> dims);

private:
	double** data_;
	double *data_2D;
};

#endif

