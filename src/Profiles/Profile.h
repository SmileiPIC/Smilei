#ifndef Profile_H
#define Profile_H


#include <vector>
#include <string>
#include <complex>
#include "SmileiMPI.h"
#include "Tools.h"
#include "PyTools.h"
#include "H5.h"
#include "Function.h"
#include "Field.h"
#include "Field3D.h"
#include "cField.h"
#include <cmath>

class Profile
{
public:
    //! Default constructor
    Profile( PyObject *, unsigned int nvariables, std::string name, Params &params, bool try_numpy=false, bool try_file=false );
    //! Cloning constructor
    Profile( Profile * );
    //! Default destructor
    ~Profile();
    
    //! Get the value of the profile at some location (spatial)
    inline double valueAt( std::vector<double> coordinates )
    {
        return function_->valueAt( coordinates );
    };
    //! Get the value of the profile at some location (temporal)
    inline double valueAt( double time )
    {
        return function_->valueAt( time );
    };
    //! Get the value of the profile at some location (spatio-temporal)
    inline double valueAt( std::vector<double> coordinates, double time )
    {
        return function_->valueAt( coordinates, time );
    };
    //! Get the complex value of the profile at some location (spatio-temporal)
    inline std::complex<double> complexValueAt( std::vector<double> coordinates, double time )
    {
        return function_->complexValueAt( coordinates, time );
    };
    
    //! Get/add the value of the profile at several locations
    //! mode = 0 : set values
    //! mode = 1 : ADD values
    //! mode = 2 : set values at given time
    //! mode = 3 : ADD values at given time
    inline void valuesAt( std::vector<Field *> &coordinates, std::vector<double> global_origin, Field &ret, int mode = 0, double time = 0 )
    {
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy_ ) {
            std::vector<PyArrayObject *> x( nvar );
            PyArrayObject *values = NULL;
            int ndim = coordinates[0]->dims().size();
            npy_intp dims[ndim];
            for( int idim=0; idim<ndim; idim++ ) {
                dims[idim] = ( npy_intp )( coordinates[0]->dims()[idim] );
            }
            // Expose arrays as numpy, and evaluate
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                x[ivar] = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( coordinates[ivar]->data() ) );
            }
            if( mode & 0b10 ) {
                values = function_->valueAt( x, time );
            } else {
                values = function_->valueAt( x );
            }
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                Py_DECREF( x[ivar] );
            }
            // Copy array to return Field3D
            double *arr = ( double * ) PyArray_GETPTR1( values, 0 );
            if( mode & 0b01 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) += arr[i];
                }
            } else {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) = arr[i];
                }
            }
            Py_DECREF( values );
        } else
#endif
        // Profile read from a file
        if( uses_file_ ) {
            std::vector<double> start( nvar );
            std::vector<double> stop ( nvar );
            std::vector<unsigned int> n = static_cast<Field3D*>(coordinates[0])->dims();
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                start[ivar]=( *coordinates[ivar] )( 0 ) - global_origin[ivar];
                stop [ivar]=( *coordinates[ivar] )( size - 1 ) - global_origin[ivar];
            }
            Field3D values = static_cast<Function_File *>( function_ )->valuesAt( start, stop, n );
            Field * v = static_cast<Field *>( &values );
            if( mode & 0b01 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) += ( *v )( i );
                }
            } else {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) = ( *v )( i );
                }
            }
        
        // Otherwise, calculate profile for each point
        } else {
            std::vector<double> x( nvar );
            if( mode == 0 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) = function_->valueAt( x );
                }
            } else if( mode == 1 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) += function_->valueAt( x );
                }
            } else if( mode == 2 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) = function_->valueAt( x, time );
                }
            } else if( mode == 3 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) += function_->valueAt( x, time );
                }
            } else {
                ERROR("valuesAt : wrong mode "<<mode);
            }
        }
    };
    
    //! Get/add the complex value of the profile at several locations
    //! mode = 0 : set values
    //! mode = 1 : ADD values
    //! mode = 2 : set values at given time
    //! mode = 3 : ADD values at given time
    inline void complexValuesAt( std::vector<Field *> &coordinates, cField &ret, int mode = 0, double time = 0 )
    {
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy_ ) {
            std::vector<PyArrayObject *> x( nvar );
            PyArrayObject *values = NULL;
            int ndim = coordinates[0]->dims().size();
            npy_intp dims[ndim];
            for( int idim=0; idim<ndim; idim++ ) {
                dims[idim] = ( npy_intp ) coordinates[0]->dims()[idim];
            }
            // Expose arrays as numpy, and evaluate
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                x[ivar] = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, ( double * )( coordinates[ivar]->data() ) );
            }
            if( mode & 0b10 ) {
                values = function_->complexValueAt( x, time );
            } else {
                values = function_->complexValueAt( x );
            }
            for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                Py_DECREF( x[ivar] );
            }
            // Copy array to return cField2D
            std::complex<double> *arr = ( std::complex<double> * ) PyArray_GETPTR1( values, 0 );
            if( mode & 0b01 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) += arr[i];
                }
            } else {
                for( unsigned int i=0; i<size; i++ ) {
                    ret( i ) = arr[i];
                }
            }
            Py_DECREF( values );
        } else
#endif
        // Profile read from a file
        if( uses_file_ ) {
            ERROR( "Profile from file not available in complex" );
        
        // Otherwise, calculate profile for each point
        } else {
            std::vector<double> x( nvar );
            if( mode == 0 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) = function_->complexValueAt( x );
                }
            } else if( mode == 1 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) += function_->complexValueAt( x );
                }
            } else if( mode == 2 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) = function_->complexValueAt( x, time );
                }
            } else if( mode == 3 ) {
                for( unsigned int i=0; i<size; i++ ) {
                    for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                        x[ivar]=( *coordinates[ivar] )( i );
                    }
                    ret( i ) += function_->complexValueAt( x, time );
                }
            } else {
                ERROR("complexValuesAt : wrong mode "<<mode);
            }
        }
    };
    
    //! Get the complex value of the profile at several locations (spatial + times)
    inline void complexValuesAtTimes( std::vector<Field *> &coordinates, Field *time, cField &ret )
    {
        unsigned int nvar = coordinates.size();
        unsigned int size = coordinates[0]->globalDims_;
#ifdef SMILEI_USE_NUMPY
        // If numpy profile, then expose coordinates as numpy before evaluating profile
        if( uses_numpy_ ) {
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
            PyArrayObject *values = function_->complexValueAt( x, t );
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
        // Profile read from a file
        if( uses_file_ ) {
            ERROR( "Profile from file not available in complex" );
        
        // Otherwise, calculate profile for each point
        } else {
            std::vector<double> x( nvar );
            double t;
            
            for( unsigned int i=0; i<size; i++ ) {
                for( unsigned int ivar=0; ivar<nvar; ivar++ ) {
                    x[ivar]=( *coordinates[ivar] )( i );
                }
                t = ( *time )( i );
                ret( i ) = function_->complexValueAt( x, t );
            }
        }
    };
    
    //! Get info on the loaded profile, to be printed later
    std::string getInfo()
    {
        std::ostringstream info( "" );
        info << nvariables_ << "D";
        
        if( ! profileName_.empty() ) {
            info << " built-in profile `" << profileName_ << "`" ;
        } else if( uses_file_ ) {
            info << " from file `" << filename_ << "`";
        } else {
            info << " user-defined function";
            if( uses_numpy_ ) {
                info << " (uses numpy)";
            }
        }
        
        if( function_ ) {
            info << function_->getInfo();
        }
        
        return info.str();
    };

private:
    
    //! Name of the profile, in the case of a built-in profile
    std::string profileName_;
    
    //! Object that holds the information on the profile function
    Function *function_;
    
    //! Number of variables for the profile function
    int nvariables_;
    
    //! Whether the profile is using numpy
    bool uses_numpy_;
    
    //! Whether the profile is taken from a file
    bool uses_file_;
    std::string filename_;
    
};//END class Profile




#endif
