
#ifndef Field3D_H
#define Field3D_H

#include "Field.h"
#include <vector>

//class Field3Dproxy {
//public:
//	friend class Field3D;
//    Field3Dproxy( double &val) : mPtr( & val ) {}
//    void operator = ( double d ) {
//    	*mPtr = d;
//    }
//    double * mPtr;
//};

class Field3D : public Field {
public:
	Field3D();
	Field3D( std::vector<unsigned int> dims );
	Field3D( std::vector<unsigned int> dims, std::string name );
	~Field3D();

	void allocateDims(std::vector<unsigned int> dims );

	inline double& operator () (unsigned int i,unsigned int j,unsigned int k) { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR("Out of limits "<< i << " " << j << " " << k));
		return data_[i][j][k]; 
	};

	inline double operator () (unsigned int i,unsigned int j,unsigned int k) const { 
		DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR("Out of limits "<< i << " " << j << " " << k));
		return data_[i][j][k]; 
	};
	

	void dump(std::vector<unsigned int> dims);

private:
	double*** data_;
	double *data_3D;
};

#endif

