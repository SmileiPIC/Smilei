#ifndef CFIELD2D_H
#define CFIELD2D_H

#include <cmath>
#include <complex>

#include <vector>

#include "Field.h"

//! class cField2D used to defined a 2d vector of Complex
class cField2D : public Field
{

public:
    //! Constructor for cField2D: no input argument
    cField2D();
    
    //! Constructor for cField2D: with the vector dimension as input argument
    cField2D( std::vector<unsigned int> dims );
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    cField2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
    //! Constructor for cField2D: with the vector dimension and filename for the dump as input argument
    cField2D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    cField2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );
    
    //! Constructor, without allocating
    cField2D( std::string name, std::vector<unsigned int> dims );
    
    //! Destructor for cField2D
    ~cField2D();
    
    //! Method used to allocate a cField2D
    void allocateDims();
    void deallocateDims();
    //! a cField2D can also be initialized win two unsigned int 
    void allocateDims(unsigned int dims1,unsigned int dims2);
    //! allocate dimensions for field2D isPrimal define if mainDim is Primal or Dual
    void allocateDims(unsigned int mainDim, bool isPrimal );
    
    inline void allocateDims(std::vector<unsigned int> dims) {
        dims_ = dims;
        allocateDims();
    };
    
    inline void allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal) {
        dims_ = dims;
        allocateDims(mainDim, isPrimal);
    };
    
    virtual void shift_x(unsigned int delta);
    
    //! Overloading of the () operator allowing to set a new value for the (i,j) element of a cField2D
    inline std::complex<double>& operator () (unsigned int i,unsigned int j) {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR(name << "Out of limits ("<< i << "," << j << ")  > (" <<dims_[0] << "," <<dims_[1] << ")" ));
        DEBUGEXEC(if (!std::isfinite(data_2D[i][j])) ERROR(name << " Not finite "<< i << "," << j << " = " << data_2D[i][j]));
        return data_2D[i][j];
    };
    
    
    //! Overloading of the () operator allowing to get the value of the (i,j) element of a cField2D
    inline std::complex<double> operator () (unsigned int i,unsigned int j) const {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR(name << "Out of limits "<< i << " " << j));
        DEBUGEXEC(if (!std::isfinite(data_2D[i][j])) ERROR(name << "Not finite "<< i << "," << j << " = " << data_2D[i][j]));
        return data_2D[i][j];
    };

    
    virtual double norm2(unsigned int istart[3][2], unsigned int bufsize[3][2]);

private:
    //! this will present the data as a 2d matrix
    std::complex<double> **data_2D;
    std::complex<double> *data_1D;
};

#endif

