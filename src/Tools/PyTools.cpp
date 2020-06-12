
#include <PyTools.h>

//! convert Python object to bool
bool PyTools::pyconvert( PyObject *py_val, bool &val )
{
    int result;
    if( py_val ) {
        result = PyObject_IsTrue( py_val );
        if( result == -1 ) {
            return false;
        }
        val = ( bool ) result;
        return true;
    }
    return false;
};

//! convert Python object to any integer type (short, uint, int, etc.)
template<typename T>
bool PyTools::pyconvert( PyObject *py_val, T &val )
{
    if( py_val && PyNumber_Check( py_val ) ) {
        PyObject *v = PyNumber_Long( py_val );
        if( PyObject_RichCompareBool( v, py_val, Py_EQ ) < 1 ) {
            return false;
        }
        val=( T ) PyLong_AsLong( v );
        Py_XDECREF( v );
        return true;
    }
    return false;
};
template bool PyTools::pyconvert<short>( PyObject *py_val, short &val );
template bool PyTools::pyconvert<unsigned int>( PyObject *py_val, unsigned int &val );
template bool PyTools::pyconvert<int>( PyObject *py_val, int &val );


//! convert Python object to double
bool PyTools::pyconvert( PyObject *py_val, double &val )
{
    if( py_val && PyNumber_Check( py_val ) ) {
        PyObject *v = PyNumber_Float( py_val );
        val=( double ) PyFloat_AsDouble( v );
        Py_XDECREF( v );
        return true;
    }
    return false;
};

//! convert Python object to complex
bool PyTools::pyconvert( PyObject *py_val, std::complex<double> &val )
{
    if( py_val && PyComplex_Check( py_val ) ) {
        val.real( PyComplex_RealAsDouble( py_val ) );
        val.imag( PyComplex_ImagAsDouble( py_val ) );
        return true;
    }
    return false;
};

//! convert Python object to string
bool PyTools::pyconvert( PyObject *py_val, std::string &val )
{
    if( py_val ) {
        PyObject *s = NULL;
        if( PyUnicode_Check( py_val ) ) {
            s = PyUnicode_AsUTF8String( py_val );    // python3
        } else if( PyBytes_Check( py_val ) ) {
            s = PyObject_Bytes( py_val );    // python2
        }
        if( s ) {
            val=std::string( PyBytes_AsString( s ) );
            Py_XDECREF( s );
            return true;
        }
    }
    return false;
};

//! check error and display message
double PyTools::get_py_result( PyObject *pyresult )
{
    checkPyError();
    double cppresult=0;
    if( pyresult ) {
        if( !py2scalar( pyresult, cppresult ) ) {
            ERROR( "A python function does not return float but " << pyresult->ob_type->tp_name );
        }
    } else {
        ERROR( "A python function raised an error" );
    }
    return cppresult;
};
std::complex<double> PyTools::get_py_result_complex( PyObject *pyresult )
{
    checkPyError();
    std::complex<double> cppresult=0;
    if( pyresult ) {
        if( !py2scalar( pyresult, cppresult ) ) {
            ERROR( "A python function does not return complex but " << pyresult->ob_type->tp_name );
        }
    } else {
        ERROR( "A python function raised an error" );
    }
    return cppresult;
};

// DECREF for vectors of python objects
void PyTools::DECREF( std::vector<PyObject *> pyvec )
{
    for( unsigned int i=0; i<pyvec.size(); i++ ) {
        if( pyvec[i] ) {
            Py_DECREF( pyvec[i] );
        }
    }
};

//! convert python list to vector of python objects
bool PyTools::py2pyvector( PyObject *py_list, std::vector<PyObject *> &py_vec )
{
    if( py_list && PySequence_Check( py_list ) ) {
        int len = PySequence_Size( py_list );
        py_vec.resize( len );
        for( int i = 0; i < len; i++ ) {
            PyObject *item = PySequence_GetItem( py_list, i );
            py_vec[i] = item;
        }
       return true; // success
    }
    checkPyError();
    return false; // error
};



//! extract a bool
void PyTools::extract( std::string name, bool &val, std::string component, int nComponent )
{
    extract( name, val, component, nComponent, "True or False" );
};

//! extract an integer-like
template< typename T>
void PyTools::extract( std::string name, T &val, std::string component, int nComponent )
{
    extract( name, val, component, nComponent, "an integer" );
};
template void PyTools::extract<short>( std::string name, short &val, std::string component, int nComponent );
template void PyTools::extract<int>( std::string name, int &val, std::string component, int nComponent );
template void PyTools::extract<unsigned int>( std::string name, unsigned int &val, std::string component, int nComponent );

//! extract a float
void PyTools::extract( std::string name, double &val, std::string component, int nComponent )
{
    extract( name, val, component, nComponent, "a float" );
};

//! extract a string
void PyTools::extract( std::string name, std::string &val, std::string component, int nComponent )
{
    extract( name, val, component, nComponent, "a string" );
};

//! extract a bool but accepts None (which returns false, and does not change the variable)
bool PyTools::extractOrNone( std::string name, bool &val, std::string component, int nComponent )
{
    return extractOrNone( name, val, component, nComponent, "True or False" );
};

//! extract an integer-like but accepts None (which returns false, and does not change the variable)
template< typename T>
bool PyTools::extractOrNone( std::string name, T &val, std::string component, int nComponent )
{
    return extractOrNone( name, val, component, nComponent, "an integer" );
};
template bool PyTools::extractOrNone<short>( std::string name, short &val, std::string component, int nComponent );
template bool PyTools::extractOrNone<int>( std::string name, int &val, std::string component, int nComponent );
template bool PyTools::extractOrNone<unsigned int>( std::string name, unsigned int &val, std::string component, int nComponent );

//! extract a float but accepts None (which returns false, and does not change the variable)
bool PyTools::extractOrNone( std::string name, double &val, std::string component, int nComponent )
{
    return extractOrNone( name, val, component, nComponent, "a float" );
};

//! extract a string but accepts None (which returns false, and does not change the variable)
bool PyTools::extractOrNone( std::string name, std::string &val, std::string component, int nComponent )
{
    return extractOrNone( name, val, component, nComponent, "a string" );
};

//! retrieve python object
PyObject * PyTools::extract_py( std::string name, std::string component, int nComponent )
{
    if( !Py_IsInitialized() ) {
        ERROR( "Python not initialized: this should not happen" );
    }
    PyObject *py_obj = PyImport_AddModule( "__main__" );
    Py_INCREF( py_obj );
    PyObject *py_return, *py_component;
    // If component requested
    if( ! component.empty() ) {
        // Get the selected component (e.g. "Species" or "Laser")
        if( PyObject_HasAttrString( py_obj, component.c_str() ) ) {
            py_component = PyObject_GetAttrString( py_obj, component.c_str() );
        } else {
            ERROR( "Component "<<component<<" is missing in the namelist" );
        }
        // If successfully found
        if( ! PySequence_Check( py_component ) ) {
            ERROR( "There is a problem with component "<<component<<". Make sure you have not overriden" );
        }
        int len = PyObject_Length( py_component );
        if( len > nComponent ) {
            py_obj = PySequence_GetItem( py_component, nComponent );
        } else {
            ERROR( "Requested " << component << " #" <<nComponent<< ", but only "<<len<<" available" );
        }
        Py_DECREF( py_component );
    }
    // Check if the variable exists and get it
    if( PyObject_HasAttrString( py_obj, name.c_str() ) ) {
        py_return = PyObject_GetAttrString( py_obj, name.c_str() );
    } else {
        ERROR( "Component "<<component<<"#"<<nComponent<<" is missing the attribute `"<<name<<"`" );
    }
    Py_DECREF( py_obj );
    return py_return;
};

//! retrieve a vector of python objects
std::vector<PyObject *> PyTools::extract_pyVec( std::string name, std::string component, int nComponent )
{
    std::vector<PyObject *> retvec;
    PyObject *py_obj = extract_py( name, component, nComponent );
    if( ! py2pyvector( py_obj, retvec ) ) {
        std::ostringstream ss( "" );
        if( component!="" ) {
            ss << "In " << component << " #" << nComponent << " ";
        }
        ERROR( ss.str() << name << " should be a list not a scalar: use [...]" );
    }
    Py_DECREF( py_obj );
    return retvec;
};

// extract 1 profile
bool PyTools::extract_pyProfile( std::string name, PyObject *&prof, std::string component, int nComponent )
{
    PyObject *py_obj = extract_py( name, component, nComponent );
    if( PyCallable_Check( py_obj ) ) {
        prof = py_obj;
        return true;
    }
    Py_XDECREF( py_obj );
    return false;
};

// extract a vector of profiles
bool PyTools::extract_pyProfiles( std::string name, std::string component, int nComponent, std::vector<PyObject *>&prof )
{
    std::vector<PyObject *> pvec = extract_pyVec( name, component, nComponent );
    for( unsigned int i=0; i<pvec.size(); i++ ) {
        if( !PyCallable_Check( pvec[i] ) ) {
            return false;
        }
    }
    prof.resize( pvec.size() );
    for( unsigned int i=0; i<pvec.size(); i++ ) {
        prof[i] = pvec[i];
    }
    return true;
};

// extract a vector of 1 or 3 profiles
bool PyTools::extract_1or3Profiles( std::string name, std::string component, int nComponent, std::vector<PyObject *>&prof )
{
    if( ! extract_pyProfiles( name, component, nComponent, prof ) ) {
        ERROR( "In "<<component<<"#"<<nComponent<<": expected a list of profiles" );
    }
    if( prof.size() == 0 ) {
        return false;
    } else if( prof.size() == 1 ) {
        prof.resize( 3 );
        prof[1] = prof[0]; Py_INCREF( prof[0] );
        prof[2] = prof[0]; Py_INCREF( prof[0] );
    } else if( prof.size() != 3 ) {
        ERROR( "In "<<component<<"#"<<nComponent<<": expected 1 or 3 profiles" );
    }
    return true;
};

// extract 2 profiles from namelist (used for laser profile)
bool PyTools::extract2Profiles( std::string varname, int ilaser, std::vector<PyObject *> &profiles )
{
    PyObject *py_obj = extract_py( varname, "Laser", ilaser );
    // Return false if None
    if( py_obj==Py_None ) {
        return false;
    }
    
    // Error if not list
    if( ! py2pyvector( py_obj, profiles ) ) {
        ERROR( "For laser #" << ilaser << ": " << varname << " must be a list of 2 profiles" );
    }
    
    // Error if wrong size
    if( profiles.size()!=2 ) {
        ERROR( "For Laser #" << ilaser << ": "<<varname<<" needs 2 profiles." );
    }
    
    // Error if not callable
    for( int i=0; i<2; i++ ) {
        if( !PyCallable_Check( profiles[i] ) ) {
            ERROR( "For Laser #" << ilaser << ": "<<varname<<"["<<i<<"] not understood" );
        }
    }
    
    return true;
};
// extract 2N profiles from namelist (used for laser profile)
bool PyTools::extract2NProfiles( std::string varname, int ilaser, std::vector<PyObject *> &profiles )
{
    PyObject *py_obj = extract_py( varname, "Laser", ilaser );
    // Return false if None
    if( py_obj==Py_None ) {
        return false;
    }
    
    // Error if not list
    if( ! py2pyvector( py_obj, profiles ) ) {
        ERROR( "For laser #" << ilaser << ": " << varname << " must be a list of 2N profiles (2 per AM)" );
    }
   

    // Error if wrong size
    if( profiles.size()%2!=0 ) {
        ERROR( "For Laser #" << ilaser << ": "<<varname<<" needs a pair number of profiles." );
    }
    
    int profile_number_of_modes = profiles.size()/2;
    // Error if not callable
    for( int i=0; i<2*profile_number_of_modes; i++ ) {
        if( !PyCallable_Check( profiles[i] ) ) {
            ERROR( "For Laser #" << ilaser << ": "<<varname<<"["<<i<<"] not understood" );
        }
    }
    
    return true;
};


//! return the number of components (see pyinit.py)
unsigned int PyTools::nComponents( std::string componentName )
{
    // Get the selected component (e.g. "Species" or "Laser")
    if( !Py_IsInitialized() ) {
        ERROR( "Python not initialized: this should not happen" );
    }
    PyObject *py_obj = PyObject_GetAttrString( PyImport_AddModule( "__main__" ), componentName.c_str() );
    checkPyError();
    Py_ssize_t retval = PyObject_Length( py_obj );
    Py_DECREF( py_obj );
    if( retval < 0 ) {
        ERROR( "Problem searching for component " << componentName );
    }
    return retval;
};


