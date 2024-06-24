#ifndef FIELD3D_H
#define FIELD3D_H

#include <cmath>

#include <vector>

#include "Field.h"
#include "Field2D.h"

class Params;
class SmileiMPI;
class Patch;

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

    //! Constructor, without allocating
    Field3D( std::string name, std::vector<unsigned int> dims );

    //! Destructor for Field3D
    ~Field3D();

    //! Method used to allocate a Field3D
    void allocateDims() override;
    void deallocateDataAndSetTo( Field* f ) override;
    //! a Field3D can also be initialized win three unsigned int
    void allocateDims( unsigned int dims1, unsigned int dims2, unsigned int dims3 );
    //! allocate dimensions for field3D isPrimal define if mainDim is Primal or Dual
    void allocateDims( unsigned int mainDim, bool isPrimal ) override;

    inline void __attribute__((always_inline)) allocateDims( std::vector<unsigned int> dims ) override
    {
        dims_ = dims;
        allocateDims();
    };

    inline void __attribute__((always_inline)) allocateDims( std::vector<unsigned int> dims, unsigned int mainDim, bool isPrimal ) override
    {
        dims_ = dims;
        allocateDims( mainDim, isPrimal );
    };

    virtual void shift_x( unsigned int delta ) override;

    //! Overloading of the () operator allowing to set a new value for the (i,j,k) element of a Field3D
    inline double & __attribute__((always_inline)) operator()( unsigned int i, unsigned int j, unsigned int k )
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k >= dims_[2] ) ERROR( name << "Out of limits & "<< i << " " << j << " " << k ) );
        return data_3D[i][j][k];
    };

    /*inline double& operator () (unsigned int i)
    {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]*dims_[2]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_3D[i])) ERROR("Not finite "<< i));
        return data_3D[i];
    };*/

    //! Overloading of the () operator allowing to get the value for the (i,j,k) element of a Field3D
    inline double __attribute__((always_inline)) operator()( unsigned int i, unsigned int j, unsigned int k ) const
    {
        DEBUGEXEC( if( i>=dims_[0] || j>=dims_[1] || k >= dims_[2] ) ERROR( name << "Out of limits "<< i << " " << j << " " << k ) );
        return data_3D[i][j][k];
    };

    /*inline double operator () (unsigned int i) const {
        DEBUGEXEC(if (i>=dims_[0]*dims_[1]*dims_[2]) ERROR("Out of limits & "<< i));
        DEBUGEXEC(if (!std::isfinite(data_3D[i])) ERROR("Not finite "<< i));
        return data_3D[i];
    };*/

    void extract_slice_yz( unsigned int ix, Field2D *field );
    void extract_slice_xz( unsigned int iy, Field2D *field );
    void extract_slice_xy( unsigned int iz, Field2D *field );

    Field* clone() override {
        auto newInstance = new Field3D(dims_);
        newInstance->name = name;
        newInstance->copyFrom(this);
        return newInstance;
    }

    //! Compute the norm2 of the field
    virtual double norm2( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override;

    //! Compute the norm2OnDevice of the field
#if defined(SMILEI_ACCELERATOR_GPU)
    virtual double norm2OnDevice( unsigned int istart[3][2], unsigned int bufsize[3][2] ) override final;
#endif

    void put( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch ) override;
    void add( Field *outField, Params &params, Patch *thisPatch, Patch *outPatch ) override;
    void get( Field  *inField, Params &params, Patch   *inPatch, Patch *thisPatch ) override;

    //double *data_3D;
    //! this will present the data as a 3d matrix
    double *** data_3D;

    void create_sub_fields  ( int iDim, int iNeighbor, int ghost_size ) override;
    void extract_fields_exch( int iDim, int iNeighbor, int ghost_size ) override;
    void inject_fields_exch ( int iDim, int iNeighbor, int ghost_size ) override;
    void extract_fields_sum ( int iDim, int iNeighbor, int ghost_size ) override;
    void inject_fields_sum  ( int iDim, int iNeighbor, int ghost_size ) override;

};

#endif
