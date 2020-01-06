#ifndef FIELD1D_H
#define FIELD1D_H

#include <cmath>
#include <vector>

#include "Field.h"
#include "Tools.h"

class Params;
class SmileiMPI;
class Patch;

//! class Field1D used to defined a 1d vector
class Field1D : public Field
{

public:
    //! Constructor for Field1D: with no input argument
    Field1D();
    
    //! Constructor for Field1D: with the vector dimension as input argument
    Field1D( std::vector<unsigned int> dims );
    
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    Field1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
    //! Constructor for Field1D: with the vector dimension and filename for the dump as input argument
    Field1D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    Field1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );
    
    //! Constructor, without allocating
    Field1D( std::string name, std::vector<unsigned int> dims );
    
    //! Destructor for Field1D
    ~Field1D();
    
    //! Method used to allocate a Field1D
    void allocateDims() override;
    void deallocateDataAndSetTo( Field* f ) override;
    //! a Field1D can also be initialized win an unsigned int
    void allocateDims( unsigned int dims1 );
    //! 1D method used to allocate Field, isPrimal define if mainDim is Primal or Dual
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
    
    //! Method to shift field in space
    void shift_x( unsigned int delta ) override;
    
    //! Overloading of the () operator allowing to set a new value for the ith element of a Field1D
    inline double &operator()( unsigned int i )
    {
        DEBUGEXEC( if( i>=dims_[0] ) ERROR( name << "Out of limits & "<< i ) );
        DEBUGEXEC( if( !std::isfinite( data_[i] ) ) ERROR( name << " not finite at i=" << i << " = " << data_[i] ) );
        return data_[i];
    };
    
    //! Overloading of the () operator allowing to get the value of the ith element of a Field1D
    inline double operator()( unsigned int i ) const
    {
        DEBUGEXEC( if( i>=dims_[0] ) ERROR( name << "Out of limits "<< i ) );
        DEBUGEXEC( if( !std::isfinite( data_[i] ) ) ERROR( name << "Not finite "<< i << " = " << data_[i] ) );
        return data_[i];
    };
    
    //! \todo What is this? (MG)
    //! \todo private/friend/modify (JD)
    //
    //! Now in Field, all arrays may be viewed as a 1D array
    //double* data_;
    
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override;
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) override;
    
private:
};


#endif
