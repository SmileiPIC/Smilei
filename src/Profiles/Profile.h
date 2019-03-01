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
#include "cField.h"
#include <cmath>

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
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy ) {
            std::vector<PyArrayObject *> x( nvar );
            int ndim = coordinates[0]->dims().size();
            npy_intp dims[ndim];
            for( int idim=0; idim<ndim; idim++ ) {
                dims[idim] = ( npy_intp )( coordinates[0]->dims()[idim] );
            }
            // Expose arrays as numpy, and evaluate
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                x[ivar] = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( coordinates[ivar]->data() ) );
            }
            PyArrayObject *values = function->valueAt( x );
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                Py_DECREF( x[ivar] );
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
            std::vector<double> x( nvar );
            for( unsigned int i=0; i<size; i++ ) {
                for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                    x[ivar]=( *coordinates[ivar] )( i );
                }
                ret( i ) = function->valueAt( x );
            }
        }
    };
    
    //! Get the value of the profile at several locations (spatial)
    inline void complexValuesAt( std::vector<Field *> &coordinates, cField &ret )
    {
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy ) {
            std::vector<PyArrayObject *> x( nvar );
            int ndim = coordinates[0]->dims().size();
            npy_intp dims[ndim];
            for( int idim=0; idim<ndim; idim++ ) {
                dims[idim] = ( npy_intp ) coordinates[0]->dims()[idim];
            }
            // Expose arrays as numpy, and evaluate
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                x[ivar] = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( coordinates[ivar]->data() ) );
            }
            PyArrayObject *values = function->complexValueAt( x );
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                Py_DECREF( x[ivar] );
            }
            // Copy array to return cField2D
            std::complex<double> *arr = ( std::complex<double> * ) PyArray_GETPTR1( values, 0 );
            for( unsigned int i=0; i<size; i++ ) {
                ret( i ) = arr[i];
            }
            Py_DECREF( values );
        } else
#endif
            // Otherwise, calculate profile for each point
        {
            std::vector<double> x( nvar );
            for( unsigned int i=0; i<size; i++ ) {
                for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                    x[ivar]=( *coordinates[ivar] )( i );
                }
                ret( i ) = function->complexValueAt( x );
            }
        }
    };
    
    
    //! Get the value of the profile at several locations (spatial)
    inline void complexValuesAt( std::vector<Field *> &coordinates, Field *time, cField &ret )
    {
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy ) {
            std::vector<PyArrayObject *> x( nvar );
            PyArrayObject *t;
            int ndim = coordinates[0]->dims().size();
            npy_intp dims[ndim];
            for( int idim=0; idim<ndim; idim++ ) {
                dims[idim] = ( npy_intp ) coordinates[0]->dims()[idim];
            }
            // Expose arrays as numpy, and evaluate
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                x[ivar] = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( coordinates[ivar]->data() ) );
            }
            t = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( time->data() ) );
            PyArrayObject *values = function->complexValueAt( x, t );
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                Py_DECREF( x[ivar] );
            }
            Py_DECREF( t );
            // Copy array to return Field3D
            std::complex<double> *arr = ( std::complex<double> * ) PyArray_GETPTR1( values, 0 );
            for( unsigned int i=0; i<size; i++ ) {
                ret( i ) = arr[i];
            }
            Py_DECREF( values );
        } else
#endif
            // Otherwise, calculate profile for each point
        {
            std::vector<double> x( nvar );
            double t;
            
            for( unsigned int i=0; i<size; i++ ) {
                for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                    x[ivar]=( *coordinates[ivar] )( i );
                }
                t = ( *time )( i );
                ret( i ) = function->complexValueAt( x, t );
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
