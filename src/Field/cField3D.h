#ifndef CFIELD3D_H
#define CFIELD3D_H

#include <cmath>
#include <complex>
#include <vector>

#include "cField.h"
#include "Patch.h"
#include "Params.h"
#include "SmileiMPI.h"

//! class cField3D used to defined a 3d vector of Complex
class cField3D : public cField
{

public:
    //! Constructor for cField3D: no input argument
    cField3D();
    
    //! Constructor for cField3D: with the vector dimension as input argument
    cField3D( std::vector<unsigned int> dims );
    //! Constructor, isPrimal define if mainDim is Primal or Dual
    cField3D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal );
    
    //! Constructor for cField3D: with the vector dimension and filename for the dump as input argument
    cField3D( std::vector<unsigned int> dims, std::string name );
    //! Constructor, isPrimal define if mainDim is Primal or Dual and a name
    cField3D( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal, std::string name );
    
    //! Constructor, without allocating
    cField3D( std::string name, std::vector<unsigned int> dims );
    
    //! Destructor for cField3D
    ~cField3D();
    
    //! Method used to allocate a cField3D
    void allocateDims() override;
    void deallocateDataAndSetTo( Field* f ) override;
    //! a cField3D can also be initialized win two unsigned int
    void allocateDims( unsigned int dims1, unsigned int dims2, unsigned int dims3 );
    //! allocate dimensions for field3D isPrimal define if mainDim is Primal or Dual
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
    
    //! Overloading of the () operator allowing to set a new value for the (i,j,k) element of a cField3D
    inline std::complex<double> &operator()( unsigned int i, unsigned int j, unsigned int k )
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k>=dims_[2] ) ERROR( name << "Out of limits ("<< i << "," << j << "," << k <<")  > (" << dims_[0] << "," << dims_[1] << "," << dims_[2] << ")" ) );
        DEBUGEXEC( if( !std::isfinite( real( data_3D[i][j][k] )+imag( data_3D[i][j][k] ) ) ) ERROR( name << " Not finite "<< i << "," << j << " = " << data_3D[i][j][k] ) );
        return data_3D[i][j][k];
    };
    
    
    //! Overloading of the () operator allowing to get the value of the (i,j,k) element of a cField3D
    inline std::complex<double> operator()( unsigned int i, unsigned int j, unsigned int k ) const
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k>=dims_[2] ) ERROR( name << "Out of limits "<< i << " " << j << " " << k ) );
        DEBUGEXEC( if( !std::isfinite( real( data_3D[i][j][k] )+imag( data_3D[i][j][k] ) ) ) ERROR( name << " Not finite "<< i << "," << j << "," << k << " = " << data_3D[i][j][k] ) );
        return data_3D[i][j][k];
    };
    
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override;
    
    inline std::complex<double> &operator()( unsigned int i )
    {
        DEBUGEXEC( if( i>=globalDims_ ) ERROR( name << " Out of limits "<< i << " < " <<dims_[0] ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[i] )+imag( cdata_[i] ) ) ) ERROR( name << " Not finite "<< i << " = " << cdata_[i] ) );
        return cdata_[i];
    };
    void dump( std::vector<unsigned int> dims ) {};
    
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) override;
    
private:
    //! this will present the data as a 2d matrix
    std::complex<double> ***data_3D;
};

#endif

