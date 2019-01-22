#ifndef Profile_H
#define Profile_H


#include <vector>
#include <string>
#include <complex>
#include "SmileiMPI.h"
#include "Tools.h"
#include "PyTools.h"
#include "Function.h"
#include "Field.h"

class Profile
{
public:
    //! Default constructor
    Profile( PyObject *, unsigned int, std::string, bool=false );
    //! Cloning constructor
    Profile( Profile * );
    //! Default destructor
    ~Profile();
    
    //! Get the value of the profile at some location (spatial)
    inline double valueAt( std::vector<double> coordinates )
    {
        return function->valueAt( coordinates );
    };
    //! Get the value of the profile at some location (temporal)
    inline double valueAt( double time )
    {
        return function->valueAt( time );
    };
    //! Get the value of the profile at some location (spatio-temporal)
    inline double valueAt( std::vector<double> coordinates, double time )
    {
        return function->valueAt( coordinates, time );
    };
    //! Get the complex value of the profile at some location (spatio-temporal)
    inline std::complex<double> complexValueAt( std::vector<double> coordinates, double time )
    {
        return function->complexValueAt( coordinates, time );
    };
    
    //! Get the value of the profile at several locations (spatial)
    inline void valuesAt( std::vector<Field *> &coordinates, Field &ret )
    {
        unsigned int ndim = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy ) {
            std::vector<PyArrayObject *> x( ndim );
            npy_intp dims[1] = {( npy_intp ) size};
            // Expose arrays as numpy, and evaluate
            for( unsigned int idim=0; idim<ndim; idim++ ) {
                x[idim] = ( PyArrayObject * )PyArray_SimpleNewFromData( 1, dims, NPY_DOUBLE, ( double * )( coordinates[idim]->data() ) );
            }
            PyArrayObject *values = function->valueAt( x );
            for( unsigned int idim=0; idim<ndim; idim++ ) {
                Py_DECREF( x[idim] );
            }
            // Copy array to return Field3D
            double *arr = ( double * ) PyArray_GETPTR1( values, 0 );
            for( unsigned int i=0; i<size; i++ ) {
                ret( i ) = arr[i];
            }
            Py_DECREF( values );
        } else
#endif
            // Otherwise, calculate profile for each point
        {
            std::vector<double> x( ndim );
            for( unsigned int i=0; i<size; i++ ) {
                // MESSAGE("  - "<<size);
                for( unsigned int idim=0; idim<ndim; idim++ ) {
                    x[idim]=( *coordinates[idim] )( i );
                }
                ret( i ) = function->valueAt( x );
            }
        }
    };
    
    //! Get info on the loaded profile, to be printed later
    inline std::string getInfo()
    {
        return info;
    };
    
    //! Name of the profile, in the case of a built-in profile
    std::string profileName;
    
private:
    //! Object that holds the information on the profile function
    Function *function;
    
    //! String containing some info on the profile
    std::string info;
    
    //! Number of variables for the profile function
    int nvariables_;
    
    //! Whether the profile is using numpy
    bool uses_numpy;
    
};//END class Profile




#endif
