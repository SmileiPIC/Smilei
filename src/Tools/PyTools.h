
// -------------------
// Some Python helper functions
// -------------------

#ifndef PYHelper_H
#define PYHelper_H


#include <Python.h>
#include <vector>
#include "Tools.h"

//! tools to query python nemlist and get back C++ values and vectors
class PyTools {
private:
    //! convert Python object to bool
    static bool pyconvert(PyObject* py_val, bool &val) {
        if (py_val && PyBool_Check(py_val)) {
            val=(py_val==Py_True);
            return true;
        }
        return false;
    }
    
    //! convert Python object to short int
    static bool pyconvert(PyObject* py_val, short int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(short int) PyInt_AsLong(py_val);
            return true;
        }
        return false;
    }
    
    //! convert Python object to unsigned int
    static bool pyconvert(PyObject* py_val, unsigned int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(unsigned int) PyInt_AsLong(py_val);
            return true;
        }
        return false;
    }
    
    //! convert Python object to int
    static bool pyconvert(PyObject* py_val, int &val) {
        if (py_val && PyInt_Check(py_val)) {
            val=(int) PyInt_AsLong(py_val);
            return true;
        }
        return false;
    }
    
    //! convert Python object to double
    static bool pyconvert(PyObject* py_val, double &val) {
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
    static bool pyconvert(PyObject* py_val, std::string &val) {
        if (py_val && PyString_Check(py_val)) {
            val=std::string(PyString_AsString(py_val));
            return true;
        }
        return false;
    }
    
    //! check error and display message
    static double get_py_result(PyObject* pyresult) {
        checkPyError();
        double cppresult=0;
        if (pyresult) {
            if(!convert(pyresult,cppresult)) {
                PyObject *ptype, *pvalue, *ptraceback;
                PyErr_Fetch(&ptype, &pvalue, &ptraceback);
                ERROR("function does not return float but " << pyresult->ob_type->tp_name);
            }
        } else {
            // we should never reach this point... something is weird
            ERROR("Function does not return a valid Python object");
        }
        return cppresult;
    }

public:

    static void openPython() {
        if (!Py_IsInitialized())
            Py_Initialize();
    }
    
    static void closePython() {
        if (Py_IsInitialized())
            Py_Finalize();
    }
    
    
    //! convert Python object to C++ value
    template <typename T>
    static bool convert(PyObject* py_vec, T &val) {
        bool retval=pyconvert(py_vec, val);
        return retval;
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
    
    //! check if there has been a python error
    static void checkPyError(bool exitOnError=false) {
        if (PyErr_Occurred()) {
            PyObject *type, *value, *traceback;
            PyErr_Fetch(&type, &value, &traceback);
            PyErr_Clear();
            
            std::string message;
            if (type) {
                type = PyObject_Str(type);
                message += PyString_AsString(type);
            }
            if (value) {
                value = PyObject_Str(value);
                message += ": ";
                message += PyString_AsString(value);
            }
            Py_XDECREF(type);
            Py_XDECREF(value);
            Py_XDECREF(traceback);
            if (exitOnError) {
                ERROR(message);
            } else {
                MESSAGE(1,"[Python] " << message);
            }
        }
    }
    
    //! run python function
    static void runPyFunction(std::string name) {
        PyObject* myFunction = PyObject_GetAttrString(PyImport_AddModule("__main__"),name.c_str());
        if (myFunction) {
            MESSAGE(1,"Calling python " << name);
            PyObject_CallFunction(myFunction,const_cast<char *>(""));
            checkPyError(true);
            Py_DECREF(myFunction);
        } else {
            MESSAGE(1,"python " << name << " function does not exists");
        }
    }

    //! get python function without arguments
    template <typename T=double>
    static T runPyFunction(std::string name) {
        T retval;
        PyObject* myFunction = PyObject_GetAttrString(PyImport_AddModule("__main__"),name.c_str());
        if (myFunction) {
            PyObject *pyresult = PyObject_CallFunction(myFunction, const_cast<char *>(""));
            retval = (T) get_py_result(pyresult);
            Py_DECREF(myFunction);
            Py_XDECREF(pyresult);
        }
        return retval;
    }
    
    //! get python function one variable    
    template <typename T=double>
    static T runPyFunction(PyObject *pyFunction, double x1) {
        PyObject *pyresult = PyObject_CallFunction(pyFunction, const_cast<char *>("d"), x1);
        T retval = (T) get_py_result(pyresult);
        Py_XDECREF(pyresult);
        return retval;
    }
    
    //! get python function two variables
    template <typename T=double>
    static T runPyFunction(PyObject *pyFunction, double x1, double x2) {
        PyObject *pyresult = PyObject_CallFunction(pyFunction, const_cast<char *>("dd"), x1, x2);
        T retval = (T) get_py_result(pyresult);
        Py_XDECREF(pyresult);
        return retval;
    }
    
    //! get python function three variables
    template <typename T=double>
    static T runPyFunction(PyObject *pyFunction, double x1, double x2, double x3) {
        PyObject *pyresult = PyObject_CallFunction(pyFunction, const_cast<char *>("ddd"), x1, x2, x3);
        T retval = (T) get_py_result(pyresult);
        Py_XDECREF(pyresult);
        return retval;
    }
    
    //! get python function four variables
    template <typename T=double>
    static T runPyFunction(PyObject *pyFunction, double x1, double x2, double x3, double x4) {
        PyObject *pyresult = PyObject_CallFunction(pyFunction, const_cast<char *>("dddd"), x1, x2, x3, x4);
        T retval = (T) get_py_result(pyresult);
        Py_XDECREF(pyresult);
        return retval;
    }
    
    //! get T from python
    template< typename T>
    static bool extract(std::string name, T &val, std::string component=std::string(""), int nComponent=0) {
        PyObject* py_val = extract_py(name,component,nComponent);
        PyTools::checkPyError();        
        return PyTools::convert(py_val,val);
    }
    
    //! extract vectors
    template< typename T>
    static bool extract(std::string name, std::vector<T> &val, std::string component=std::string(""), int nComponent=0) {
        std::vector<PyObject*> py_val = extract_pyVec(name,component,nComponent);
        if (py_val.size())
            return PyTools::convert(py_val,val);
        return false;
    }
    
    //! retrieve python object
    static PyObject* extract_py(std::string name, std::string component=std::string(""), int nComponent=0) {
        if (name.find(" ")!= std::string::npos || component.find(" ")!= std::string::npos) {
            WARNING("asking for [" << name << "] [" << component << "] : it has whitespace inside: please fix the code");
        }
        if (!Py_IsInitialized()) ERROR("Python not initialized: this should not happend");
        PyObject *py_obj=PyImport_AddModule("__main__");
        // If component requested
        if (!component.empty()) {
            // Get the selected component (e.g. "Species" or "Laser")
            py_obj = PyObject_GetAttrString(py_obj,component.c_str());
            PyTools::checkPyError();
            // Error if not found
            if (!py_obj) ERROR("Component "<<component<<" not found in namelist");
            // If successfully found
            int len = PyObject_Length(py_obj);
            if (len > nComponent) {
                py_obj = PySequence_GetItem(py_obj, nComponent);
            } else {
                ERROR("Requested " << component << " #" <<nComponent<< ", but only "<<len<<" available");
            }
        }
        PyObject *py_return=PyObject_GetAttrString(py_obj,name.c_str());
        PyTools::checkPyError();
        return py_return;
        
    }
    
    static void toProfile(PyObject*& myPy) {
        double val;
        // If the profile is only a double, then convert to a constant function
        if( PyTools::convert(myPy, val) ) {
            // Extract the function "constant"
            PyObject* constantFunction = PyTools::extract_py("constant");
            // Create the argument which has the value of the profile
            PyObject* arg = PyTuple_New(1);
            PyTuple_SET_ITEM(arg, 0, PyFloat_FromDouble(val));
            // Create the constant anonymous function
            myPy = PyObject_Call(constantFunction, arg, NULL);
            PyTools::checkPyError();
        }
    }

    static bool extract_pyProfile(std::string name, PyObject*& myPy, std::string component=std::string(""), int nComponent=0) {
        PyObject* myPytmp=extract_py(name,component,nComponent);
        toProfile(myPytmp);
        if (PyCallable_Check(myPytmp)) {
            myPy=myPytmp;
            return true;
        }
        return false;
    }
    
    //! retrieve a vector of python objects
    static std::vector<PyObject*> extract_pyVec(std::string name, std::string component=std::string(""), int nComponent=0) {
        std::vector<PyObject*> retvec;
        PyObject* py_obj = extract_py(name,component,nComponent);
        if (py_obj) {
            if (!PyTuple_Check(py_obj) && !PyList_Check(py_obj)) {
                retvec.push_back(py_obj);
                WARNING(name << " should be a list or tuple, not a scalar : fix it");
            } else {
                PyObject* seq = PySequence_Fast(py_obj, "expected a sequence");
                int len = PySequence_Size(py_obj);
                retvec.resize(len);
                for (int i = 0; i < len; i++) {
                    PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
                    retvec[i]=item;
                }
                Py_DECREF(seq);
            }
        }
        PyTools::checkPyError();
        return retvec;
    }
    
    // extract 3 profiles from namelist (used for part mean velocity and temperature)
    static void extract3Profiles(std::string varname, int ispec, PyObject*& profx, PyObject*& profy, PyObject*& profz )
    {
        std::vector<PyObject*> pvec = PyTools::extract_pyVec(varname,"Species",ispec);
        for (unsigned int i=0;i<pvec.size();i++) {
            PyTools::toProfile(pvec[i]);
        }
        if ( pvec.size()==1 ) {
            profx =  profy =  profz = pvec[0];
        } else if (pvec.size()==3) {
            profx = pvec[0];
            profy = pvec[1];
            profz = pvec[2];
        } else {
            ERROR("For species #" << ispec << ", "<<varname<<" needs 1 or 3 components.");
        }
    }

    //! return the number of components (see pyinit.py)
    static int nComponents(std::string componentName) {
        // Get the selected component (e.g. "Species" or "Laser")
        if (!Py_IsInitialized()) ERROR("Python not initialized: this should not happend");
        PyObject *py_obj = PyObject_GetAttrString(PyImport_AddModule("__main__"),componentName.c_str());
        PyTools::checkPyError();
        int retval = PyObject_Length(py_obj);
        return retval;
    }
    
    
};

#endif
