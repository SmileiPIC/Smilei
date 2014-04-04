#ifndef FIELD1D_H
#define FIELD1D_H

#include <cmath>
#include <vector>

#include "Field.h"
#include "Tools.h"

//! class Field1D used to defined a 1d vector
class Field1D : public Field
{

public:
    //! Constructor for Field1D: with no input argument
    Field1D();

    //! Constructor for Field1D: with the vector dimension as input argument
    Field1D( std::vector<unsigned int> dims );
    Field1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );

    //! Constructor for Field1D: with the vector dimension and filename for the dump as input argument
    Field1D( std::vector<unsigned int> dims, std::string name );
    Field1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );

    //! Destructor for Field1D
    ~Field1D();

    //! Method used to allocate a Field1D
    void allocateDims(std::vector<unsigned int> dims);
    void allocateDims(std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal);

    //! Method used to dump the data contained in a Field1D
    void dump(std::vector<unsigned int> dims);

    //! Overloading of the () operator allowing to set a new value for the ith element of a Field1D
    inline double& operator () (unsigned int i)
    {
        DEBUGEXEC(if (i>=dims_[0]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_[i])) ERROR("Not finite "<< i << " = " << data_[i]));
        return data_[i];
    };

    //! Overloading of the () operator allowing to get the value of the ith element of a Field1D
    inline double operator () (unsigned int i) const
    {
        DEBUGEXEC(if (i>=dims_[0]) ERROR("Out of limits "<< i));
        DEBUGEXEC(if (!std::isfinite(data_[i])) ERROR("Not finite "<< i << " = " << data_[i]));
        return data_[i];
    };

    //! \todo What is this? (MG)
    //! \todo private/friend/modify SmileiMPI* (JD)
    // 
    //! Now in Field, all arrays may be viewed as a 1D array
    //double* data_;

private:
};


#endif
