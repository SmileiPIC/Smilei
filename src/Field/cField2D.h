#ifndef CFIELD2D_H
#define CFIELD2D_H

#include <cmath>
#include <complex>

#include <vector>

#include "cField.h"

//! class cField2D used to defined a 2d vector of Complex
class cField2D : public cField
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
    void allocateDims() override;
    void deallocateDataAndSetTo( Field* f ) override;
    //! a cField2D can also be initialized win two unsigned int
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
    
    //! Overloading of the () operator allowing to set a new value for the (i,j) element of a cField2D
    inline std::complex<double> &operator()( unsigned int i, unsigned int j )
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] ) ERROR( name << "Out of limits ("<< i << "," << j << ")  > (" <<dims_[0] << "," <<dims_[1] << ")" ) );
        DEBUGEXEC(if ( !std::isfinite( real(data_2D[i][j])+imag(data_2D[i][j]) ) ) ERROR(name << " Not finite "<< i << "," << j << " = " << data_2D[i][j] ));
        //DEBUGEXEC( if( std::abs( data_2D[i][j] ) > 1.2 ) ERROR( name << " Greater than 1.2 "<< i << "," << j << " = " << data_2D[i][j] ) );
        return data_2D[i][j];
    };
    
    
    //! Overloading of the () operator allowing to get the value of the (i,j) element of a cField2D
    inline std::complex<double> operator()( unsigned int i, unsigned int j ) const
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] ) ERROR( name << "Out of limits "<< i << " " << j ) );
        DEBUGEXEC(if (!std::isfinite(real(data_2D[i][j])+imag(data_2D[i][j]))) ERROR(name << " Not finite "<< i << "," << j << " = " << data_2D[i][j] ));
        //DEBUGEXEC( if( std::abs( data_2D[i][j] ) > 1.2 ) ERROR( name << " Greater than 1.2 "<< i << "," << j << " = " << data_2D[i][j] ) );
        return data_2D[i][j];
    };
    
    
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override;
    
    inline std::complex<double> &operator()( unsigned int i )
    {
        DEBUGEXEC( if( i>=globalDims_ ) ERROR( name << " Out of limits "<< i << " < " <<dims_[0] ) );
        DEBUGEXEC( if( !std::isfinite( real( cdata_[i] )+imag( cdata_[i] ) ) ) ERROR( name << " Not finite "<< i << " = " << cdata_[i] ) );
        return cdata_[i];
    };
    
    //! method used to put all entry of a field at a given value val
    void put_to( double val ) override
    {
        if( cdata_ )
            for( unsigned int i=0; i<globalDims_; i++ ) {
                cdata_[i] = val;
            }
    }
    void copyFrom( Field *from_field )
    {
        cField2D *from_cfield = static_cast<cField2D *>( from_field );
        DEBUGEXEC( if( globalDims_!=from_field->globalDims_ ) ERROR( "Field size do not match "<< name << " " << from_field->name ) );
        for( unsigned int i=0; i< globalDims_; i++ ) {
            ( *this )( i )=( *from_cfield )( i );
        }
    }
    
    void put( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, SmileiMPI *smpi, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, SmileiMPI *smpi, Patch   *inPatch, Patch *thisPatch ) override;
    
    //! this will present the data as a 2d matrix
    std::complex<double> **data_2D;

private:

};

#endif

