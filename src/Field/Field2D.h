#ifndef FIELD2D_H
#define FIELD2D_H

#include <cmath>

#include <vector>

#include "Field.h"

class Params;
class SmileiMPI;
class Patch;

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
    
    //! Constructor, without allocating
    Field2D( std::string name, std::vector<unsigned int> dims );
    
    //! Destructor for Field2D
    ~Field2D();
    
    //! Method used to allocate a Field2D
    void allocateDims() override;
    void deallocateDataAndSetTo( Field* f ) override;
    //! a Field2D can also be initialized win two unsigned int
    void allocateDims( unsigned int dims1, unsigned int dims2 );
    //! allocate dimensions for field2D isPrimal define if mainDim is Primal or Dual
    void allocateDims( unsigned int mainDim, bool isPrimal ) override;
    
    inline void allocateDims( std::vector<unsigned int> dims ) override
    {
        dims_ = dims;
        allocateDims();
    };
    
    inline void allocateDims( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) override
    {
        dims_ = dims;
        allocateDims( mainDim, isPrimal );
    };
    
    virtual void shift_x( unsigned int delta ) override;
    
    //! Overloading of the () operator allowing to set a new value for the (i,j) element of a Field2D
    inline double &operator()( unsigned int i, unsigned int j )
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] ) ERROR( name << "Out of limits ("<< i << "," << j << ")  > (" <<dims_[0] << "," <<dims_[1] << ")" ) );
        DEBUGEXEC( if( !std::isfinite( data_2D[i][j] ) ) ERROR( name << " Not finite "<< i << "," << j << " = " << data_2D[i][j] ) );
        return data_2D[i][j];
    };
    
    /*inline double& operator () (unsigned int i) {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_2D[i])) ERROR("Not finite "<< i));
        return data_2D[i];
    };*/
    
    //! Overloading of the () operator allowing to get the value of the (i,j) element of a Field2D
    inline double operator()( unsigned int i, unsigned int j ) const
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] ) ERROR( name << "Out of limits "<< i << " " << j ) );
        DEBUGEXEC( if( !std::isfinite( data_2D[i][j] ) ) ERROR( name << "Not finite "<< i << "," << j << " = " << data_2D[i][j] ) );
        return data_2D[i][j];
    };
    
    /*inline double operator () (unsigned int i) const {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_2D[i])) ERROR("Not finite "<< i));
        return data_2D[i];
    };*/
    
    //double** data_;
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override;
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) override;
    
    //!\todo{Comment what are these stuffs (MG for JD)}
    //double *data_2D;
    //! this will present the data as a 2d matrix
    double **data_2D;
    
};

#endif

