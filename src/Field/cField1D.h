#ifndef CFIELD1D_H
#define CFIELD1D_H

#include <cmath>
#include <complex>

#include <vector>

#include "cField.h"
#include "Tools.h"

class Params;
class SmileiMPI;
class Patch;

//! class Field1D used to defined a 1d vector
class cField1D : public cField
{

public:
    //! Constructor for cField1D: with no input argument
    cField1D();
    
    //! Constructor for cField1D: with the vector dimension as input argument
    cField1D( std::vector<unsigned int> dims );
    
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    cField1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
    //! Constructor for cField1D: with the vector dimension and filename for the dump as input argument
    cField1D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    cField1D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );
    
    //! Constructor, without allocating
    cField1D( std::string name, std::vector<unsigned int> dims );
    
    //! Destructor for Field1D
    ~cField1D();
    
    //! Method used to allocate a Field1D
    void allocateDims();
    void deallocateDataAndSetTo( Field* f ) override;
    //! a Field1D can also be initialized win an unsigned int
    void allocateDims( unsigned int dims1 );
    //! 1D method used to allocate Field, isPrimal define if mainDim is Primal or Dual
    void allocateDims( unsigned int mainDim, bool isPrimal );
    
    inline void allocateDims( std::vector<unsigned int> dims )
    {
        dims_ = dims;
        allocateDims();
    };
    
    inline void allocateDims( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal )
    {
        dims_ = dims;
        allocateDims( mainDim, isPrimal );
    };
    
    //! Method to shift field in space
    void shift_x( unsigned int delta );
    
    //! Overloading of the () operator allowing to set a new value for the ith element of a Field1D
    inline std::complex<double> &operator()( unsigned int i )
    {
        DEBUGEXEC( if( i>=dims_[0] ) ERROR( name << "Out of limits & "<< i ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[i] )+imag( cdata_[i] ) ) ) ERROR( name << " not finite at i=" << i << " = " << cdata_[i] ) );
        return cdata_[i];
    };
    
    //! Overloading of the () operator allowing to get the value of the ith element of a Field1D
    inline std::complex<double> operator()( unsigned int i ) const
    {
        DEBUGEXEC( if( i>=dims_[0] ) ERROR( name << "Out of limits "<< i ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[i] )+imag( cdata_[i] ) ) ) ERROR( name << "Not finite "<< i << " = " << cdata_[i] ) );
        return cdata_[i];
    };
    
    //! \todo What is this? (MG)
    //! \todo private/friend/modify (JD)
    //
    //! Now in Field, all arrays may be viewed as a 1D array
    //double* data_;
    
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] );
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) override;
    
private:
};


#endif
