
// -------------------
// Some Python helper functions
// -------------------

#ifndef PYHelper_H
#define PYHelper_H


#include <Python.h>
#include "Tools.h"

//! tools to convert python values to C++ values and vectors
class PyTools {    
public:

    //! convert Python object to bool
    static bool convert(PyObject* py_val, bool &val) {
        if (py_val && PyBool_Check(py_val)) {
            val=(py_val==Py_True);
            return true;
        }
        return false;        
    }
    
    //! convert Python object to short int
    static bool convert(PyObject* py_val, short int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(short int) PyInt_AsLong(py_val);
            return true;
        }
        return false;        
    }
    
    //! convert Python object to unsigned int
    static bool convert(PyObject* py_val, unsigned int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(unsigned int) PyInt_AsLong(py_val);
            return true;
        }
        return false;        
    }
    
    //! convert Python object to int
    static bool convert(PyObject* py_val, int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(int) PyInt_AsLong(py_val);
            return true;
        }
        return false;        
    }
    
    //! convert Python object to double
    static bool convert(PyObject* py_val, double &val) {
        if(py_val) {
            if (PyFloat_Check(py_val)) {
                val = PyFloat_AsDouble(py_val);
                return true;
            } else if (PyInt_Check(py_val)) {
                val=(double) PyInt_AsLong(py_val);
                return true;
            }
        }
        return false;        
    }

    //! convert Python object to string
    static bool convert(PyObject* py_val, std::string &val) {
        if (py_val && PyString_Check(py_val)) {
            val=std::string(PyString_AsString(py_val));
            return true;
        }
        return false;
    }
    
    //! convert vector of Python objects to vector of C++ values
    template <typename T>
    static bool convert(std::vector<PyObject*> py_vec, std::vector<T> &val) {
        bool retval=true;
        val.resize(py_vec.size());
        for (unsigned int i=0;i<py_vec.size();i++) {
            bool thisval=convert(py_vec[i], val[i]);
            if (thisval==false) retval=false;
        }
        return retval;
    }
    
};

#endif