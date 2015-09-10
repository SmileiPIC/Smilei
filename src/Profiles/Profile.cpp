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
