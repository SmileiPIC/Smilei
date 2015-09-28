#include "ElectroMagn.h"
#include "Profile.h"

using namespace std;

// Preliminary functions
// that evaluate a python function with various numbers of arguments
double Evaluate1var(PyObject * fun, std::vector<double> x_cell) {
    return PyTools::runPyFunction(fun, x_cell[0]);
}
double Evaluate2var(PyObject * fun, std::vector<double> x_cell) {
    return PyTools::runPyFunction(fun, x_cell[0], x_cell[1]);
}
double Evaluate3var(PyObject * fun, std::vector<double> x_cell) {
    return PyTools::runPyFunction(fun, x_cell[0], x_cell[1], x_cell[2]);
}


// Default constructor.
Profile::Profile(PyObject* pp, unsigned int nvariables) :
py_profile(pp) {
    
    if (!PyCallable_Check(py_profile)) {
        ERROR("Profile can't be called");
    }
    
    PyObject* repr = PyObject_Repr(py_profile);
    DEBUG(string(PyString_AsString(repr)));
    Py_XDECREF(repr);

    
    repr=PyObject_Str(py_profile);
    DEBUG(string(PyString_AsString(repr)));
    Py_XDECREF(repr);
    
    PyObject* inspect=PyImport_ImportModule("inspect");

    PyTools::checkPyError();
    PyObject *tuple = PyObject_CallMethod(inspect,const_cast<char *>("getargspec"),const_cast<char *>("(O)"),py_profile);

    PyObject *arglist = PyTuple_GetItem(tuple,0);
    int size = PyObject_Size(arglist);
    if (size != (int)nvariables) {
        WARNING ("Profile takes "<< size <<" variables but it is crated with " << nvariables);
        string args("");
        for (int i=0; i<size; i++){
            PyObject *arg=PyList_GetItem(arglist,i);
            
            PyObject* repr = PyObject_Repr(arg);
            args+=string(PyString_AsString(repr))+" ";
            Py_XDECREF(repr);
        }
        MESSAGE("Profile vars: "<<args);
    }
    
    
    Py_XDECREF(tuple);
    Py_XDECREF(inspect);
    

    
    if      ( nvariables == 1 ) Evaluate = &Evaluate1var;
    else if ( nvariables == 2 ) Evaluate = &Evaluate2var;
    else if ( nvariables == 3 ) Evaluate = &Evaluate3var;
    else {
        ERROR("A profile has been defined with unsupported number of variables");
    }
}


double Profile::valueAt(vector<double> x_cell) {
    return (*Evaluate)(py_profile, x_cell);
}
