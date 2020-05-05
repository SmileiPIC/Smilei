#ifndef CFIELD_H
#define CFIELD_H

#include <cmath>
#include <complex>
#include <vector>

#include "Field.h"
#include "Params.h"
#include "SmileiMPI.h"

class Patch;

//! class cField used to defined an array of Complex
class cField : public Field
{

public:
    //! Constructor for cField: no input argument
    cField() : Field()
    {
    };
    //! Constructor for cField: with the vector dimension as input argument
    cField( std::vector<unsigned int> dims ) : Field( dims )
    {
    };
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    cField( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) : Field( dims, mainDim, isPrimal )
    {
    };
    //! Constructor for cField: with the vector dimension and filename for the dump as input argument
    cField( std::vector<unsigned int> dims, std::string name_in ) : Field( dims, name_in )
    {
    };
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    cField( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name_in ) : Field( dims, mainDim, isPrimal, name_in )
    {
    };
    
    //! Destructor for cField
    virtual ~cField()
    {
    };
    
    //! Method used to allocate a cField
    virtual void allocateDims() = 0;
    virtual void deallocateDataAndSetTo( Field* f ) = 0;
    //! a cField can also be initialized win two unsigned int
//    void allocateDims(unsigned int dims1,unsigned int dims2,unsigned int dims3);
//    //! allocate dimensions for field3D isPrimal define if mainDim is Primal or Dual
//    void allocateDims(unsigned int mainDim, bool isPrimal );

    virtual void allocateDims( std::vector<unsigned int> dims ) = 0;
    virtual void allocateDims( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) = 0;
    
    virtual void shift_x( unsigned int delta ) = 0;
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) = 0;
    
    inline std::complex<double> &operator()( unsigned int i )
    {
        DEBUGEXEC( if( i>=globalDims_ ) ERROR( name << " Out of limits "<< i << " < " <<dims_[0] ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[i] )+imag( cdata_[i] ) ) ) ERROR( name << " Not finite "<< i << " = " << cdata_[i] ) );
        return cdata_[i];
    };

    //! 2D reference access to the linearized array (with check in DEBUG mode)
    inline std::complex<double> &operator()( unsigned int i, unsigned int j )
    {
        int unsigned idx = i*dims_[1]+j;
        DEBUGEXEC( if( idx>=globalDims_ ) ERROR( "Out of limits & "<< i << " " << j ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[idx] )+imag( cdata_[idx] ) ) ) ERROR( "Not finite "<< i << " " << j << " = " << cdata_[idx] ) );
        return cdata_[idx];
    };
    //! 2D access to the linearized array (with check in DEBUG mode)
    inline std::complex<double> operator()( unsigned int i, unsigned int j ) const
    {
        unsigned int idx = i*dims_[1]+j;
        DEBUGEXEC( if( idx>=globalDims_ ) ERROR( "Out of limits "<< i << " " << j ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[idx] )+imag( cdata_[idx] ) ) ) ERROR( "Not finite "<< i << " " << j << " = " << cdata_[idx] ) );
        return cdata_[idx];
    };
    
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) = 0;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) = 0;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) = 0;
    
    std::complex<double> *cdata_;
protected:
};

#endif

