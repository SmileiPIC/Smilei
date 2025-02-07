
// -------------------
// Some Python helper functions
// -------------------

#ifndef PYHelper_H
#define PYHelper_H


#if defined(_POSIX_C_SOURCE) || defined(_XOPEN_SOURCE)
//#warning check for fix here: https://bugs.python.org/issue17120
#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif
#endif

#include <Python.h>
#include <typeinfo>
#include <vector>
#include <sstream>
#include <complex>
#include "Tools.h"

#ifdef SMILEI_USE_NUMPY
#define PY_ARRAY_UNIQUE_SYMBOL SMILEI_ARRAY_API
#ifndef SMILEI_IMPORT_ARRAY
#define NO_IMPORT_ARRAY
#endif
#include <numpy/arrayobject.h>
#endif

// The following definitions are required for some interaction between python and openmp
// Specifically, it is needed for an openmp loop that includes calls to python
// Usage:
//  SMILEI_PY_SAVE_MASTER_THREAD
//  #pragma omp for
//  for(....) {
//      SMILEI_PY_ACQUIRE_GIL
//      python calls
//      SMILEI_PY_RELEASE_GIL
//  }
//  SMILEI_PY_RESTORE_MASTER_THREAD
#if PY_MAJOR_VERSION > 2 && PY_MINOR_VERSION > 11 
    #define SMILEI_PY_SAVE_MASTER_THREAD  PyThreadState *_save = NULL; _Pragma("omp master") if( Py_IsInitialized() ){ _save = PyEval_SaveThread();}  _Pragma("omp barrier")
    #define SMILEI_PY_RESTORE_MASTER_THREAD _Pragma("omp master") if( Py_IsInitialized() ){ PyEval_RestoreThread( _save );}  _Pragma("omp barrier")
    #define SMILEI_PY_ACQUIRE_GIL PyGILState_STATE _state = PyGILState_Ensure();
    #define SMILEI_PY_RELEASE_GIL PyGILState_Release( _state );
#else
    #define SMILEI_PY_SAVE_MASTER_THREAD
    #define SMILEI_PY_RESTORE_MASTER_THREAD
    #define SMILEI_PY_ACQUIRE_GIL _Pragma("omp critical") {
    #define SMILEI_PY_RELEASE_GIL }
#endif


//! tools to query python nemlist and get back C++ values and vectors
class PyTools
{
private:
    //! convert Python object to bool
    static bool pyconvert( PyObject *py_val, bool &val );

    //! convert Python object to any integer type (short, uint, int, etc.)
    template<typename T>
    static bool pyconvert( PyObject *py_val, T &val );

    //! convert Python object to double
    static bool pyconvert( PyObject *py_val, double &val );

    //! convert Python object to complex
    static bool pyconvert( PyObject *py_val, std::complex<double> &val );

    //! convert Python object to string
    static bool pyconvert( PyObject *py_val, std::string &val );

    //! check error and display message
    template <typename T>
    static T get_py_result( PyObject *pyresult )
    {
        checkPyError();
        T cppresult = 0;
        if( pyresult ) {
            if( !py2scalar( pyresult, cppresult ) ) {
                std::string type_name;
                double a_double;
                std::complex<double> a_complex;
                if( typeid(T) == typeid(a_double) ) {
                    type_name = "float";
                } else if( typeid(T) == typeid(a_complex) ) {
                    type_name = "complex";
                } else {
                    type_name = "???";
                }
                ERROR( "A python function does not return "<< type_name <<" but " << pyresult->ob_type->tp_name );
            }
        } else {
            ERROR( "A python function raised an error" );
        }
        return cppresult;
    }

    // DECREF for vectors of python objects
    static void DECREF( std::vector<PyObject *> pyvec );

public:

    static void openPython()
    {
        if( !Py_IsInitialized() ) {
            Py_Initialize();
        }
    }

    static void closePython()
    {
        if( Py_IsInitialized() ) {
            Py_Finalize();
        }
    }

    static std::string python_version()
    {
        std::string version;
        PyObject *platform = PyImport_ImportModule( "platform" );
        PyObject *pyversion = PyObject_CallMethod( platform, const_cast<char *>( "python_version" ), NULL );
        pyconvert( pyversion, version );
        Py_XDECREF( pyversion );
        Py_XDECREF( platform );
        PyErr_Print();
        return version;
    }

    //! convert Python object to C++ value
    template <typename T>
    static bool py2scalar( PyObject *py_val, T &val )
    {
        bool retval=pyconvert( py_val, val );
        return retval;
    }

    //! convert vector of Python objects to vector of C++ values
    template <typename T>
    static bool pyvector2vector( std::vector<PyObject *> py_vec, std::vector<T> &val )
    {
        bool retval = true;
        val.resize( py_vec.size() );
        for( unsigned int i=0; i<py_vec.size(); i++ ) {
            if( ! py2scalar( py_vec[i], val[i] ) ) {
                retval = false;
            }
        }
        return retval;
    }

    //! convert python list to vector of python objects
    static bool py2pyvector( PyObject *py_list, std::vector<PyObject *> &py_vec );

    //! convert python list to vector of c++ values
    template <typename T>
    static bool py2vector( PyObject *py_list, std::vector<T> &vec )
    {
        std::vector<PyObject *> py_vec;
        bool ret = false;
        if( py2pyvector( py_list, py_vec ) ) {
            if( pyvector2vector( py_vec, vec ) ) {
                ret = true;
            }
            DECREF( py_vec );
        }
        return ret;
    };

    //! check if there has been a python error
    static bool checkPyError( bool exitOnError=false, bool print=true )
    {
        if( PyErr_Occurred() ) {
            if( exitOnError || print ) {
                PyObject *type, *value, *traceback;
                PyErr_Fetch( &type, &value, &traceback );

                std::ostringstream message( "" );
                if( type ) {
                    std::string str( "" );
                    PyObject *tn = PyObject_GetAttrString( type, "__name__" );
                    pyconvert( tn, str );
                    Py_XDECREF( tn );
                    message << str;
                }
                if( value ) {
                    message << ": " << repr( value );
                }
                Py_XDECREF( type );
                Py_XDECREF( value );
                Py_XDECREF( traceback );
                if( exitOnError ) {
                    ERROR( message.str() );
                } else if( print ) {
                    int rk(0);
                    MPI_Comm_rank( MPI_COMM_WORLD, &rk );
                    std::ostringstream t("");
                    t << "  On rank "<< rk <<" [Python] " << message.str() << std::endl;
                    std::cout << t.str();
                }
            }
            PyErr_Clear();
            return true;
        }
        return false;
    }

    //! run void python function
    static void runPyFunction( std::string name )
    {
        PyObject *myFunction = PyObject_GetAttrString( PyImport_AddModule( "__main__" ), name.c_str() );
        if( myFunction ) {
            MESSAGE( 1, "Calling python " << name );
            PyObject_CallFunction( myFunction, const_cast<char *>( "" ) );
            checkPyError( true );
            Py_DECREF( myFunction );
        } else {
            MESSAGE( 1, "python " << name << " function does not exist" );
        }
    }

    //! run typed python function without arguments
    template <typename T=double>
    static T runPyFunction( std::string name, std::string component=std::string( "" ) )
    {
        T retval( 0 );
        PyObject *py_obj=PyImport_AddModule( "__main__" );
        if( !component.empty() ) {
            py_obj = PyObject_GetAttrString( py_obj, component.c_str() );
            checkPyError();
        }
        PyObject *myFunction = PyObject_GetAttrString( py_obj, name.c_str() );
        if( myFunction ) {
            PyObject *pyresult = PyObject_CallFunction( myFunction, const_cast<char *>( "" ) );
            retval = get_py_result<T>( pyresult );
            Py_DECREF( myFunction );
            Py_XDECREF( pyresult );
        }
        return retval;
    }

    //! run typed python function with one argument
    template <typename T=double>
    static T runPyFunction( PyObject *pyFunction, double x1 )
    {
        PyObject *pyresult = PyObject_CallFunction( pyFunction, const_cast<char *>( "d" ), x1 );
        T retval = get_py_result<T>( pyresult );
        Py_XDECREF( pyresult );
        return retval;
    }

    //! run typed python function with two arguments
    template <typename T=double>
    static T runPyFunction( PyObject *pyFunction, double x1, double x2 )
    {
        PyObject *pyresult = PyObject_CallFunction( pyFunction, const_cast<char *>( "dd" ), x1, x2 );
        T retval = get_py_result<T>( pyresult );
        Py_XDECREF( pyresult );
        return retval;
    }

    //! run typed python function with three arguments
    template <typename T=double>
    static T runPyFunction( PyObject *pyFunction, double x1, double x2, double x3 )
    {
        PyObject *pyresult = PyObject_CallFunction( pyFunction, const_cast<char *>( "ddd" ), x1, x2, x3 );
        T retval = get_py_result<T>( pyresult );
        Py_XDECREF( pyresult );
        return retval;
    }

    //! run typed python function with four arguments
    template <typename T=double>
    static T runPyFunction( PyObject *pyFunction, double x1, double x2, double x3, double x4 )
    {
        PyObject *pyresult = PyObject_CallFunction( pyFunction, const_cast<char *>( "dddd" ), x1, x2, x3, x4 );
        T retval = get_py_result<T>( pyresult );
        Py_XDECREF( pyresult );
        return retval;
    }

    //! extract a bool
    static void extract( std::string name, bool &val, std::string component=std::string( "" ), int nComponent=0 );
    //! extract an integer-like
    template< typename T>
    static void extract( std::string name, T &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract a float
    static void extract( std::string name, double &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract a string
    static void extract( std::string name, std::string &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract any scalar
    template< typename T>
    static void extract( std::string name, T &val, std::string component, int nComponent, std::string testMessage )
    {
        PyObject *py_val = extract_py( name, component, nComponent );
        checkPyError();
        if( ! py2scalar( py_val, val ) ) {
            std::string print_type;
            py2scalar( PyObject_Str( PyObject_Type( py_val ) ), print_type );
            ERROR( "In "<<component<<"#"<<nComponent<<": `"<<name<<"` should be "<<testMessage<<" but is "<<print_type );
        }
    }
    
    //! extract a bool but accepts None (which returns false, and does not change the variable)
    static bool extractOrNone( std::string name, bool &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract an integer-like but accepts None (which returns false, and does not change the variable)
    template< typename T>
    static bool extractOrNone( std::string name, T &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract a float but accepts None (which returns false, and does not change the variable)
    static bool extractOrNone( std::string name, double &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract a string but accepts None (which returns false, and does not change the variable)
    static bool extractOrNone( std::string name, std::string &val, std::string component=std::string( "" ), int nComponent=0 );

    //! extract any scalar but accepts None (which returns false, and does not change the variable)
    template< typename T>
    static bool extractOrNone( std::string name, T &val, std::string component, int nComponent, std::string testMessage )
    {
        PyObject *py_val = extract_py( name, component, nComponent );
        checkPyError();
        if( py_val == Py_None ) {
            Py_XDECREF( py_val );
            return false;
        } else if( ! py2scalar( py_val, val ) ) {
            std::string print_type;
            py2scalar( PyObject_Str( PyObject_Type( py_val ) ), print_type );
            ERROR( "In "<<component<<"#"<<nComponent<<": `"<<name<<"` should be "<<testMessage<<" but is "<<print_type );
        }
        return true;
    }
    
    //! extract vector
    template< typename T>
    static bool extractV( std::string name, std::vector<T> &val, std::string component, int nComponent=0 )
    {
        std::vector<PyObject *> py_val = extract_pyVec( name, component, nComponent );
        if( py_val.size() ) {
            return pyvector2vector( py_val, val );
        }
        return false;
    }
    
    //! extract vector of vectors
    template< typename T>
    static bool extractVV( std::string name, std::vector<std::vector<T> > &val, std::string component, int nComponent=0 )
    {
        std::vector<PyObject *> py_val = extract_pyVec( name, component, nComponent );
        if( py_val.size() == 0 ) {
            return false;
        }
        val.resize( 0 );
        for( unsigned int i=0; i<py_val.size(); i++ ) {
            std::vector<T> vec;
            if( ! py2vector( py_val[i], vec ) ) {
                ERROR( name << " should be a list of lists" );
            }
            val.push_back( vec );
        }
        return true;
    }
    
    //! retrieve python object
    static PyObject *extract_py( std::string name, std::string component=std::string( "" ), int nComponent=0 );

    //! retrieve a vector of python objects
    static std::vector<PyObject *> extract_pyVec( std::string name, std::string component=std::string( "" ), int nComponent=0 );

    // extract 1 profile
    static bool extract_pyProfile( std::string name, PyObject *&prof, std::string component=std::string( "" ), int nComponent=0 );

    // extract a vector of profiles
    static bool extract_pyProfiles( std::string name, std::string component, int nComponent, std::vector<PyObject *>&prof );

    // extract a vector of 1 or N profiles
    static bool extract_1orNProfiles( unsigned int numberOfProfiles, std::string name, std::string component, int nComponent, std::vector<PyObject *>&prof );

    // extract 2 profiles from namelist (used for laser profile)
    static bool extract2Profiles( std::string varname, int ilaser, std::vector<PyObject *> &profiles );

    // extract 2N profiles from namelist (used for laser profile)
    static bool extract2NProfiles( std::string varname, int ilaser, std::vector<PyObject *> &profiles );


    //! return the number of components (see pyinit.py)
    static unsigned int nComponents( std::string componentName );

    //! Get an object's attribute ( int, double, etc.)
    template <typename T>
    static bool getAttr( PyObject *object, std::string attr_name, T &value )
    {
        bool success = false;
        if( PyObject_HasAttrString( object, attr_name.c_str() ) ) {
            PyObject *py_value = PyObject_GetAttrString( object, attr_name.c_str() );
            success = py2scalar( py_value, value );
            Py_XDECREF( py_value );
        }
        return success;
    }
    
    //! Get an object's attribute for lists
    template <typename T>
    static bool getAttr( PyObject *object, std::string attr_name, std::vector<T> &vec )
    {
        bool success = false;
        if( PyObject_HasAttrString( object, attr_name.c_str() ) ) {
            PyObject *py_list = PyObject_GetAttrString( object, attr_name.c_str() );
            success = py2vector( py_list, vec );
            Py_XDECREF( py_list );
        }
        return success;
    }
    
    //! Get an object's attribute for lists of list
    template <typename T>
    static bool getAttr( PyObject *object, std::string attr_name, std::vector<std::vector<T> > &vec )
    {
        if( PyObject_HasAttrString( object, attr_name.c_str() ) ) {
            // Get the list of lists
            PyObject *py_list = PyObject_GetAttrString( object, attr_name.c_str() );
            // Convert to vector of lists
            std::vector<PyObject *> py_vec;
            if( !py2pyvector( py_list, py_vec ) ) {
                Py_XDECREF( py_list );
                return false;
            }
            Py_XDECREF( py_list );
            // For each list, convert to vector
            std::vector<T> v;
            std::vector<std::vector<T> > vv( 0 );
            for( unsigned int i=0; i<py_vec.size(); i++ ) {
                if( ! py2vector( py_vec[i], v ) ) {
                    DECREF( py_vec );
                    return false;
                }
                vv.push_back( v );
            }
            DECREF( py_vec );
            vec = vv;
            return true;
        }
        //This should never happen
        return false;
    }
    
    //! Get an object's repr
    static std::string repr( PyObject *obj )
    {
        PyObject *s = PyObject_Str( obj );
        std::string str( "" );
        pyconvert( s, str );
        Py_XDECREF( s );
        return str;
    }

    // Set a python variable to itime so that it can be accessed at run-time
    inline static void setIteration( int itime )
    {
        PyObject *Main = PyObject_GetAttrString( PyImport_AddModule( "__main__" ), "Main" );
        PyObject *iteration = PyLong_FromLong( itime );
        PyObject_SetAttrString( Main, "iteration", iteration );
        Py_DECREF( iteration );
        Py_DECREF( Main );
    }

    // Test whether a python object is a function
    // Return the number of arguments
    // If there are python *args, return -1
    // If failure, return -2
    static int function_nargs( PyObject *obj )
    {
        if( ! PyCallable_Check( obj ) ) {
            return -2;
        }
        PyObject *inspect = PyImport_ImportModule( "inspect" );
        checkPyError();
#if PY_MAJOR_VERSION < 3
        PyObject *tuple = PyObject_CallMethod( inspect, const_cast<char *>( "getargspec" ), const_cast<char *>( "(O)" ), obj );
#else
        PyObject *tuple = PyObject_CallMethod( inspect, const_cast<char *>( "getfullargspec" ), const_cast<char *>( "(O)" ), obj );
#endif
        PyObject *arglist = PyTuple_GetItem( tuple, 0 ); // list of function arguments
        PyObject *vararg = PyTuple_GetItem( tuple, 1 ); // name of *args
        int nargs = PyObject_Size( arglist );
        bool varargs = vararg != Py_None;
        Py_XDECREF( tuple );
        Py_XDECREF( inspect );
        if( varargs ) {
            return -1;
        } else {
            return nargs;
        }
    }

    //! Finds if there is a species with a given name
    static bool isSpecies( std::string name );
};

#endif
