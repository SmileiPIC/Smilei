
#ifndef FIELD1D_H
#define FIELD1D_H

#include "Field.h"
#include <vector>
#include <cmath>

#include "Tools.h"

class Field1D : public Field {
public:
	Field1D();
	Field1D( std::vector<unsigned int> dims );
	Field1D( std::vector<unsigned int> dims, std::string name );
	~Field1D();
	
	inline double& operator () (unsigned int i) {
		DEBUGEXEC(if (i>=dims_[0]) ERROR("Out of limits "<< i));
		DEBUGEXEC(if (!std::isfinite(data_[i])) ERROR("Not finite "<< i << " = " << data_[i]));
		return data_[i];
	};

	void allocateDims(std::vector<unsigned int> dims);

	inline double  operator () (unsigned int i) const { 
		DEBUGEXEC(if (i>=dims_[0]) ERROR("Out of limits "<< i));
		DEBUGEXEC(if (!std::isfinite(data_[i])) ERROR("Not finite "<< i << " = " << data_[i]));
		return data_[i]; 
	};

	void dump(std::vector<unsigned int> dims);
private:
	double* data_;
};

#endif

