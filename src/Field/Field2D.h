#ifndef FIELD2D_H
#define FIELD2D_H

#include <cmath>

#include <vector>

#include "Field.h"

//! class Field2D used to defined a 2d vector
class Field2D : public Field
{

public:
    //! Constructor for Field2D: no input argument
    Field2D();

    //! Constructor for Field2D: with the vector dimension as input argument
    Field2D( std::vector<unsigned int> dims );
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );

    //! Constructor for Field2D: with the vector dimension and filename for the dump as input argument
    Field2D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    Field2D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );

    //! Destructor for Field2D
    ~Field2D();

    //! Method used to allocate a Field2D
    void allocateDims(std::vector<unsigned int> dims );
    //! a Field2D can also be initialized win two unsigned int 
    void allocateDims(unsigned int dims1,unsigned int dims2);
    //! allocate dimensions for field2D isPrimal define if mainDim is Primal or Dual
    void allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );

    //! Method used to dump the data contained in a Field2D
    void dump(std::vector<unsigned int> dims);
    virtual void shift_x(unsigned int delta);

    //! Overloading of the () operator allowing to set a new value for the (i,j) element of a Field2D
    inline double& operator () (unsigned int i,unsigned int j) {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR(name << "Out of limits & "<< i << " " << j));
        DEBUGEXEC(if (!std::isfinite(data_2D[i][j])) ERROR(name << " Not finite "<< i << "," << j << " = " << data_2D[i][j]));
        return data_2D[i][j];
    };

    /*inline double& operator () (unsigned int i) {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_2D[i])) ERROR("Not finite "<< i));
        return data_2D[i];
    };*/

    //! Overloading of the () operator allowing to get the value of the (i,j) element of a Field2D
    inline double operator () (unsigned int i,unsigned int j) const {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1]) ERROR(name << "Out of limits "<< i << " " << j));
        DEBUGEXEC(if (!std::isfinite(data_2D[i][j])) ERROR(name << "Not finite "<< i << "," << j << " = " << data_2D[i][j]));
        return data_2D[i][j];
    };

    /*inline double operator () (unsigned int i) const {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_2D[i])) ERROR("Not finite "<< i));
        return data_2D[i];
    };*/

    //double** data_;
    //! this will present the data as a 2d matrix
    double **data_2D;

    //virtual double computeNRJ(unsigned int shift, unsigned int** istart, unsigned int** bufsize) {return 0.;};
    virtual double computeNRJ(unsigned int shift, unsigned int istart[3][2], unsigned int bufsize[3][2]);

private:
    //!\todo{Comment what are these stuffs (MG for JD)}
    //double *data_2D;
};

#endif

