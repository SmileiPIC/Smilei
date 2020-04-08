#include <cmath>

#include "Profile.h"
#include "PyTools.h"

using namespace std;



// Default constructor.
Profile::Profile( PyObject *py_profile, unsigned int nvariables, string name, bool try_numpy ) :
    profileName( "" ),
    nvariables_( nvariables ),
    uses_numpy( false )
{
    ostringstream info_( "" );
    info_ << nvariables_ << "D";
    
    if( !PyCallable_Check( py_profile ) ) {
        ERROR( "Profile `"<<name<<"`: not a function" );
    }
    
    // In case the function was created in "pyprofiles.py", then we transform it
    //  in a "hard-coded" function
    if( PyObject_HasAttrString( py_profile, "profileName" ) ) {
    
        PyTools::getAttr( py_profile, "profileName", profileName );
        
        info_ << " built-in profile `" << profileName << "`" ;
        
        if( profileName == "constant" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Constant1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Constant2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Constant3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: constant() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "trapezoidal" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Trapezoidal1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Trapezoidal2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Trapezoidal3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: trapezoidal() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "gaussian" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Gaussian1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Gaussian2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Gaussian3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: gaussian() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "polygonal" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Polygonal1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Polygonal2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Polygonal3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: polygonal() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "cosine" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Cosine1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Cosine2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Cosine3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: cosine() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "polynomial" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_Polynomial1D( py_profile );
            } else if( nvariables_ == 2 ) {
                function = new Function_Polynomial2D( py_profile );
            } else if( nvariables_ == 3 ) {
                function = new Function_Polynomial3D( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: polynomial() profile defined only in 1D, 2D or 3D" );
            }
            
        } else if( profileName == "tconstant" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimeConstant( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tconstant() profile is only for time" );
            }
            
        } else if( profileName == "ttrapezoidal" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimeTrapezoidal( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: ttrapezoidal() profile is only for time" );
            }
            
        } else if( profileName == "tgaussian" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimeGaussian( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tgaussian() profile is only for time" );
            }
            
        } else if( profileName == "tpolygonal" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimePolygonal( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tpolygonal() profile is only for time" );
            }
            
        } else if( profileName == "tcosine" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimeCosine( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tcosine() profile is only for time" );
            }
            
        } else if( profileName == "tpolynomial" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimePolynomial( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tpolynomial() profile is only for time" );
            }
            
        } else if( profileName == "tsin2plateau" ) {
        
            if( nvariables_ == 1 ) {
                function = new Function_TimeSin2Plateau( py_profile );
            } else {
                ERROR( "Profile `"<<name<<"`: tsin2plateau() profile is only for time" );
            }
        } else {
            
            ERROR( "Undefined profile "<<profileName );
        
        }
        
    }
    
    // Otherwise (if the python profile cannot be hard-coded) ....
    else {
        string message;
        
#ifdef  __DEBUG
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
            WARNING( "Profile `" << name << "` takes "<< nargs <<" arguments but requires " << nvariables_ );
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
                PyObject *ret( nullptr );
                if( nvariables_ == 1 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, NULL );
                } else if( nvariables_ == 2 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, a, NULL );
                } else if( nvariables_ == 3 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, a, a, NULL );
                } else if( nvariables_ == 4 ) {
                    ret = PyObject_CallFunctionObjArgs( py_profile, a, a, a, a, NULL );
                }
#ifdef  __DEBUG
                DEBUG( "Profile `"<<name<<"`: try numpy array of dimension " << ndim );
                PyTools::checkPyError( false, true );
#else
                PyTools::checkPyError( false, false );
#endif
                Py_DECREF( a );
                if( ret
                        && PyArray_Check( ret ) // must be a numpy array
                        && PyArray_ISNUMBER( ( PyArrayObject * )ret ) // must be an array of floats
                        && PyArray_SIZE( ( PyArrayObject * )ret ) == numel ) { // must have the same size as arguments
                    uses_numpy = true;
                }
                if( ret ) {
                    Py_DECREF( ret );
                }
                if( uses_numpy ) {
                    break;
                }
            }
            if( uses_numpy ) {
                DEBUG( "Profile `"<<name<<"`: accepts numpy arrays of dimension " << ndim );
            } else {
                DEBUG( "Profile `"<<name<<"`: does not seem to accept numpy arrays (and will be slow)" );
            }
        }
#endif
        if( !uses_numpy ) {
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
            function = new Function_Python1D( py_profile );
        } else if( nvariables_ == 2 ) {
            function = new Function_Python2D( py_profile );
        } else if( nvariables_ == 3 ) {
            function = new Function_Python3D( py_profile );
        } else if( nvariables_ == 4 ) {
            function = new Function_Python4D( py_profile );
        }
        
        info_ << " user-defined function";
        if( try_numpy ) {
            if( uses_numpy ) {
                info_ << " (uses numpy)";
            } else {
                info_ << " (does not use numpy)";
            }
        }
    }
    info = info_.str();
    if (function) {
        info += function->getInfo();
    }
}


// Cloning constructor
Profile::Profile( Profile *p )
{
    profileName = p->profileName;
    nvariables_ = p->nvariables_;
    info        = p->info       ;
    uses_numpy  = p->uses_numpy ;
    if( profileName != "" ) {
        if( profileName == "constant" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Constant1D( static_cast<Function_Constant1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Constant2D( static_cast<Function_Constant2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Constant3D( static_cast<Function_Constant3D *>( p->function ) );
            }
        } else if( profileName == "trapezoidal" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Trapezoidal1D( static_cast<Function_Trapezoidal1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Trapezoidal2D( static_cast<Function_Trapezoidal2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Trapezoidal3D( static_cast<Function_Trapezoidal3D *>( p->function ) );
            }
        } else if( profileName == "gaussian" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Gaussian1D( static_cast<Function_Gaussian1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Gaussian2D( static_cast<Function_Gaussian2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Gaussian3D( static_cast<Function_Gaussian3D *>( p->function ) );
            }
        } else if( profileName == "polygonal" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Polygonal1D( static_cast<Function_Polygonal1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Polygonal2D( static_cast<Function_Polygonal2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Polygonal3D( static_cast<Function_Polygonal3D *>( p->function ) );
            }
        } else if( profileName == "cosine" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Cosine1D( static_cast<Function_Cosine1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Cosine2D( static_cast<Function_Cosine2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Cosine3D( static_cast<Function_Cosine3D *>( p->function ) );
            }
        } else if( profileName == "polynomial" ) {
            if( nvariables_ == 1 ) {
                function = new Function_Polynomial1D( static_cast<Function_Polynomial1D *>( p->function ) );
            } else if( nvariables_ == 2 ) {
                function = new Function_Polynomial2D( static_cast<Function_Polynomial2D *>( p->function ) );
            } else if( nvariables_ == 3 ) {
                function = new Function_Polynomial3D( static_cast<Function_Polynomial3D *>( p->function ) );
            }
        } else if( profileName == "tconstant" ) {
            function = new Function_TimeConstant( static_cast<Function_TimeConstant *>( p->function ) );
        } else if( profileName == "ttrapezoidal" ) {
            function = new Function_TimeTrapezoidal( static_cast<Function_TimeTrapezoidal *>( p->function ) );
        } else if( profileName == "tgaussian" ) {
            function = new Function_TimeGaussian( static_cast<Function_TimeGaussian *>( p->function ) );
        } else if( profileName == "tpolygonal" ) {
            function = new Function_TimePolygonal( static_cast<Function_TimePolygonal *>( p->function ) );
        } else if( profileName == "tcosine" ) {
            function = new Function_TimeCosine( static_cast<Function_TimeCosine *>( p->function ) );
        } else if( profileName == "tpolynomial" ) {
            function = new Function_TimePolynomial( static_cast<Function_TimePolynomial *>( p->function ) );
        } else if( profileName == "tsin2plateau" ) {
            function = new Function_TimeSin2Plateau( static_cast<Function_TimeSin2Plateau *>( p->function ) );
        }
    } else {
        if( nvariables_ == 1 ) {
            function = new Function_Python1D( static_cast<Function_Python1D *>( p->function ) );
        } else if( nvariables_ == 2 ) {
            function = new Function_Python2D( static_cast<Function_Python2D *>( p->function ) );
        } else if( nvariables_ == 3 ) {
            function = new Function_Python3D( static_cast<Function_Python3D *>( p->function ) );
        } else if( nvariables_ == 4 ) {
            function = new Function_Python4D( static_cast<Function_Python4D *>( p->function ) );
        }
    }
}



Profile::~Profile()
{
    delete function;
}
