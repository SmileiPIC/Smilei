#include <cmath>

#include "Profile.h"
#include "PyTools.h"

using namespace std;



// Default constructor.
Profile::Profile( PyObject *py_profile, unsigned int nvariables, string name, Params &params, bool try_numpy, bool try_file, bool time_variable ) :
    profileName_( "" ),
    nvariables_( nvariables ),
    uses_numpy_( false ),
    uses_file_( false ),
    filename_( "" )
{
    // In case the function was created in "pyprofiles.py", then we transform it
    //  in a "hard-coded" function
    if( PyObject_HasAttrString( py_profile, "profileName" ) ) {
    
        PyTools::getAttr( py_profile, "profileName", profileName_ );
        
        if( profileName_ == "constant" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Constant1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Constant2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Constant3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: constant() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "trapezoidal" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Trapezoidal1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Trapezoidal2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Trapezoidal3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: trapezoidal() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "gaussian" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Gaussian1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Gaussian2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Gaussian3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: gaussian() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "polygonal" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Polygonal1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Polygonal2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Polygonal3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: polygonal() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "cosine" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Cosine1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Cosine2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Cosine3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: cosine() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "polynomial" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_Polynomial1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Polynomial2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Polynomial3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: polynomial() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName_ == "tconstant" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimeConstant( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tconstant() profile is only for time" );
            }
            
        } else if( profileName_ == "ttrapezoidal" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimeTrapezoidal( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: ttrapezoidal() profile is only for time" );
            }
            
        } else if( profileName_ == "tgaussian" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimeGaussian( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tgaussian() profile is only for time" );
            }
            
        } else if( profileName_ == "tpolygonal" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimePolygonal( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tpolygonal() profile is only for time" );
            }
            
        } else if( profileName_ == "tcosine" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimeCosine( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tcosine() profile is only for time" );
            }
            
        } else if( profileName_ == "tpolynomial" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimePolynomial( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tpolynomial() profile is only for time" );
            }
            
        } else if( profileName_ == "tsin2plateau" ) {
        
            if( nvariables_ == 1 ) {
                function_ = new Function_TimeSin2Plateau( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tsin2plateau() profile is only for time" );
            }
        } else {
            
            ERROR( "Undefined profile "<<profileName_ );
        
        }
        
    }
    
    // If the profile is a python function (not hard-coded)
    else if( PyCallable_Check( py_profile ) ) {
        
#ifdef  __DEBUG
        string message;
        // Check how the profile looks like
        PyObject *repr = PyObject_Repr( py_profile );
        PyTools::py2scalar( repr, message );
        MESSAGE( message );
        Py_XDECREF( repr );
        
        repr = PyObject_Str( py_profile );
        PyTools::py2scalar( repr, message );
        MESSAGE( message );
        Py_XDECREF( repr );
#endif
        
        // Compatibility with python 2.7
        // Some class functions (for SmileiSingleton) may not be static.
        // In that case, we look for the __func__ attribute.
        if( PyObject_HasAttrString( py_profile, "__func__" ) ) {
            PyObject *func = PyObject_GetAttrString( py_profile, "__func__" );
            Py_DECREF( py_profile );
            py_profile = func;
        }
        
        // Verify that the profile has the right number of arguments
        int nargs = PyTools::function_nargs( py_profile );
        if( nargs == -2 ) {
            ERROR( "Profile `" << name << "` does not seem to be callable" );
        }
        if( nargs >= 0  && nargs != ( int ) nvariables_ ) {
            WARNING( "Profile `" << name << "` has "<< nargs <<" arguments but requires " << nvariables_ );
        }
        if( nvariables_<1 || nvariables_>4 ) {
            ERROR( "Profile `"<<name<<"`: defined with unsupported number of variables (" << nvariables_ << ")" );
        }
        
        
        // Verify that the profile transforms a float in a float
#ifdef SMILEI_USE_NUMPY
        if( try_numpy ) {
            // If numpy available, verify that the profile accepts numpy arguments
            // We test 2 options : the arrays dimension equal to nvariables or nvariables-1
            unsigned int ndim;
            for( ndim=nvariables-1; ndim<=nvariables; ndim++ ) {
                unsigned int numel = pow( 2, ndim );
                double test_value[numel];
                for( unsigned int i=0; i<numel; i++ ) {
                    test_value[i] = 0.;
                }
                npy_intp dims[ndim];
                for( unsigned int i=0; i<ndim; i++ ) {
                    dims[i] = 2;
                }
                PyArrayObject *a = ( PyArrayObject * )PyArray_SimpleNewFromData( ndim, dims, NPY_DOUBLE, &test_value );
                PyObject *t = time_variable ? PyFloat_FromDouble(0.) : (PyObject *)a;
                PyObject *ret( nullptr );
                if( nvariables_ == 1 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, NULL );
                } else if( nvariables_ == 2 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, t, NULL );
                } else if( nvariables_ == 3 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, a, t, NULL );
                } else if( nvariables_ == 4 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, a, a, t, NULL );
                }
#ifdef  __DEBUG
                DEBUG( "Profile `"<<name<<"`: try numpy array of dimension " << ndim );
                PyTools::checkPyError( false, true );
#else
                PyTools::checkPyError( false, false );
#endif
                if( time_variable ) Py_DECREF( t );
                Py_DECREF( a );
                if( ret && PyArray_Check( ret ) // must be a numpy array
                        && PyArray_ISNUMBER( ( PyArrayObject * )ret ) // must be an array of floats
                        && PyArray_SIZE( ( PyArrayObject * )ret ) == numel ) { // must have the same size as arguments
                    uses_numpy_ = true;
                }
                if( ret ) {
                    Py_DECREF( ret );
                }
                if( uses_numpy_ ) {
                    break;
                }
            }
            if( uses_numpy_ ) {
                DEBUG( "Profile `"<<name<<"`: accepts numpy arrays of dimension " << ndim );
            } else {
                DEBUG( "Profile `"<<name<<"`: does not seem to accept numpy arrays (and will be slow)" );
            }
        }
#endif
        if( !uses_numpy_ ) {
            // Otherwise, try a float
            PyObject *z = PyFloat_FromDouble( 0. );
            PyObject *ret( nullptr );
            if( nvariables_ == 1 ) {
                ret = PyObject_CallFunctionObjArgs( py_profile, z, NULL );
            } else if( nvariables_ == 2 ) {
                ret = PyObject_CallFunctionObjArgs( py_profile, z, z, NULL );
            } else if( nvariables_ == 3 ) {
                ret = PyObject_CallFunctionObjArgs( py_profile, z, z, z, NULL );
            } else if( nvariables_ == 4 ) {
                ret = PyObject_CallFunctionObjArgs( py_profile, z, z, z, z, NULL );
            }
            PyTools::checkPyError();
            Py_DECREF( z );
            if( !ret || !PyNumber_Check( ret ) ) {
                ERROR( "Profile `"<<name<<"`: does not seem to return a correct value" );
            }
            if( ret ) {
                Py_DECREF( ret );
            }
        }
        
        // Assign the evaluating function, which depends on the number of arguments
        if( nvariables_ == 1 ) {
            function_ = new Function_Python1D( py_profile );
        } else if( nvariables_  == 2 ) {
            function_ = new Function_Python2D( py_profile );
        } else if( nvariables_  == 3 ) {
            function_ = new Function_Python3D( py_profile );
        } else if( nvariables_  == 4 ) {
            function_ = new Function_Python4D( py_profile );
        }
    
    // If the profile is a string (hdf5 file)
    } else if( PyTools::py2scalar( py_profile, filename_ ) ) {
        
        if( ! try_file ) {
            ERROR( "Profile `"<<name<<"`: cannot be from a file" );
        }
        
        // The filename_ should formed like myfile.h5/path/to/dataset
        // We separate the last slash (/) to isolate the dataset name
        size_t i = filename_.find_last_of( "/" );
        if( i == std::string::npos || i == filename_.size() - 1 ) {
            ERROR( "Profile `"<<name<<"`: could not find a dataset name in the path" );
        }
        string dataset_name = filename_.substr( i + 1 );
        string path = filename_.substr( 0, i );
        
        // Open the file + group
        H5Read * file = new H5Read( path );
        
        // Create the function
        function_ = new Function_File( path, dataset_name, file, params.cell_length );
        uses_file_ = true;
        
    } else {
        ERROR( "Profile `"<<name<<"`: should be a function or a string" );
    }
}


// Cloning constructor
Profile::Profile( Profile *p )
{
    profileName_ = p->profileName_;
    nvariables_ = p->nvariables_;
    uses_numpy_  = p->uses_numpy_ ;
    uses_file_ = p->uses_file_;
    filename_ = p->filename_;
    
    if( profileName_ != "" ) {
        if( profileName_ == "constant" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Constant1D( static_cast<Function_Constant1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Constant2D( static_cast<Function_Constant2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Constant3D( static_cast<Function_Constant3D *>( p->function_ ) );
            }
        } else if( profileName_ == "trapezoidal" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Trapezoidal1D( static_cast<Function_Trapezoidal1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Trapezoidal2D( static_cast<Function_Trapezoidal2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Trapezoidal3D( static_cast<Function_Trapezoidal3D *>( p->function_ ) );
            }
        } else if( profileName_ == "gaussian" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Gaussian1D( static_cast<Function_Gaussian1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Gaussian2D( static_cast<Function_Gaussian2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Gaussian3D( static_cast<Function_Gaussian3D *>( p->function_ ) );
            }
        } else if( profileName_ == "polygonal" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Polygonal1D( static_cast<Function_Polygonal1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Polygonal2D( static_cast<Function_Polygonal2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Polygonal3D( static_cast<Function_Polygonal3D *>( p->function_ ) );
            }
        } else if( profileName_ == "cosine" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Cosine1D( static_cast<Function_Cosine1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Cosine2D( static_cast<Function_Cosine2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Cosine3D( static_cast<Function_Cosine3D *>( p->function_ ) );
            }
        } else if( profileName_ == "polynomial" ) {
            if( nvariables_ == 1 ) {
                function_ = new Function_Polynomial1D( static_cast<Function_Polynomial1D *>( p->function_ ) );
            } else if( nvariables_ == 2 ) {
                function_ = new Function_Polynomial2D( static_cast<Function_Polynomial2D *>( p->function_ ) );
            } else if( nvariables_ == 3 ) {
                function_ = new Function_Polynomial3D( static_cast<Function_Polynomial3D *>( p->function_ ) );
            }
        } else if( profileName_ == "tconstant" ) {
            function_ = new Function_TimeConstant( static_cast<Function_TimeConstant *>( p->function_ ) );
        } else if( profileName_ == "ttrapezoidal" ) {
            function_ = new Function_TimeTrapezoidal( static_cast<Function_TimeTrapezoidal *>( p->function_ ) );
        } else if( profileName_ == "tgaussian" ) {
            function_ = new Function_TimeGaussian( static_cast<Function_TimeGaussian *>( p->function_ ) );
        } else if( profileName_ == "tpolygonal" ) {
            function_ = new Function_TimePolygonal( static_cast<Function_TimePolygonal *>( p->function_ ) );
        } else if( profileName_ == "tcosine" ) {
            function_ = new Function_TimeCosine( static_cast<Function_TimeCosine *>( p->function_ ) );
        } else if( profileName_ == "tpolynomial" ) {
            function_ = new Function_TimePolynomial( static_cast<Function_TimePolynomial *>( p->function_ ) );
        } else if( profileName_ == "tsin2plateau" ) {
            function_ = new Function_TimeSin2Plateau( static_cast<Function_TimeSin2Plateau *>( p->function_ ) );
        }
    } else if( uses_file_ ) {
        function_ = new Function_File( static_cast<Function_File *>( p->function_ ) );
    } else {
        if( nvariables_ == 1 ) {
            function_ = new Function_Python1D( static_cast<Function_Python1D *>( p->function_ ) );
        } else if( nvariables_ == 2 ) {
            function_ = new Function_Python2D( static_cast<Function_Python2D *>( p->function_ ) );
        } else if( nvariables_ == 3 ) {
            function_ = new Function_Python3D( static_cast<Function_Python3D *>( p->function_ ) );
        } else if( nvariables_ == 4 ) {
            function_ = new Function_Python4D( static_cast<Function_Python4D *>( p->function_ ) );
        }
    }
}



Profile::~Profile()
{
    delete function_;
}

//! Get/add the value of the profile at several locations
//! mode = 0 : set values
//! mode = 1 : ADD values
//! mode = 2 : set values at given time
//! mode = 3 : ADD values at given time
void Profile::valuesAt( std::vector<Field *> &coordinates, std::vector<double> global_origin, Field &ret, int mode, double time )
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
}

//! Get/add the complex value of the profile at several locations
//! mode = 0 : set values
//! mode = 1 : ADD values
//! mode = 2 : set values at given time
//! mode = 3 : ADD values at given time
void Profile::complexValuesAt( std::vector<Field *> &coordinates, cField &ret, int mode, double time )
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
}

//! Get the complex value of the profile at several locations (spatial + times)
void Profile::complexValuesAtTimes( std::vector<Field *> &coordinates, Field *time, cField &ret )
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
}