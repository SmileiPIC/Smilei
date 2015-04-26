#ifndef FIELD3D_H
#define FIELD3D_H

#include <cmath>

#include <vector>

#include "Field.h"

//! class Field3D used to defined a 3d vector
class Field3D : public Field
{

public:
    //! Constructor for Field3D: no input argument
    Field3D();

    //! Constructor for Field2D: with the vector dimension as input argument
    Field3D( std::vector<unsigned int> dims );
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    Field3D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );

    //! Constructor for Field2D: with the vector dimension and filename for the dump as input argument
    Field3D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    Field3D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );

    //! Destructor for Field3D
    ~Field3D();

    //! Method used to allocate a Field3D
    void allocateDims(std::vector<unsigned int> dims );
    //! a Field3D can also be initialized win three unsigned int 
	void allocateDims(unsigned int dims1,unsigned int dims2,unsigned int dims3);
    //! allocate dimensions for field3D isPrimal define if mainDim is Primal or Dual
	void allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );

    //! Method used to dump the data contained in a Field3D
    void dump(std::vector<unsigned int> dims);
    virtual void shift_x(unsigned int delta);

    //! Overloading of the () operator allowing to set a new value for the (i,j,k) element of a Field3D
    inline double& operator () (unsigned int i,unsigned int j,unsigned int k)
    {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR(name << "Out of limits & "<< i << " " << j << " " << k));
        return data_3D[i][j][k];
    };

    /*inline double& operator () (unsigned int i)
    {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]*dims_[2]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_3D[i])) ERROR("Not finite "<< i));
        return data_3D[i];
    };*/

    //! Overloading of the () operator allowing to get the value for the (i,j,k) element of a Field3D
    inline double operator () (unsigned int i,unsigned int j,unsigned int k) const {
        DEBUGEXEC(if (i>=dims_[0] || j>=dims_[1] || k >= dims_[2]) ERROR(name << "Out of limits "<< i << " " << j << " " << k));
        return data_3D[i][j][k];
    };

    /*inline double operator () (unsigned int i) const {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]*dims_[2]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_3D[i])) ERROR("Not finite "<< i));
        return data_3D[i];
    };*/

    //virtual double computeNRJ(unsigned int shift, unsigned int** istart, unsigned int** bufsize) {return 0.;};
    virtual double computeNRJ(unsigned int shift, unsigned int istart[3][2], unsigned int bufsize[3][2]);

private:
    //!\todo{Comment what are these stuffs (MG for JD)}
    //double*** data_;
    //double *data_3D;
    
    //! this will present the data as a 3d matrix
    double ***data_3D;
};

#endif
