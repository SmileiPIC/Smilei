
// -------------------
// Some Python helper functions
// -------------------

#ifndef PYHelper_H
#define PYHelper_H


#include <Python.h>
#include "PicParams.h"
#include "Tools.h"

class PyHelper {
    
public:
    
    static double py_eval_profile(ProfileStructure prof, double val0) {
        std::vector<double> val(1);
        val[0]=val0;
        return py_eval_profile(prof,val);
    }
    
    static double py_eval_profile(ProfileStructure prof, double val0, double val1) {
        std::vector<double> val(2);
        val[0]=val0;
        val[1]=val1;
        return py_eval_profile(prof,val);
    }

    static double py_eval_profile(ProfileStructure prof, double val0, double val1, double val2) {
        std::vector<double> val(3);
        val[0]=val0;
        val[1]=val1;
        val[2]=val2;
        return py_eval_profile(prof,val);
    }
    
    static double py_eval_profile(ProfileStructure prof, double val0, double val1, double val2, double val3) {
        std::vector<double> val(4);
        val[0]=val0;
        val[1]=val1;
        val[2]=val2;
        val[3]=val3;
        return py_eval_profile(prof,val);
    }
    
    static double py_eval_profile(ProfileStructure prof, std::vector<double> val) {
        PyObject *myvals=PyTuple_New(val.size()+PyTuple_Size(prof.py_args));
        for (int i=0;i<val.size();i++)
            PyTuple_SetItem(myvals, i, PyFloat_FromDouble(val[i]));
        for (int i=0;i<PyTuple_Size(prof.py_args);i++)
            PyTuple_SetItem(myvals, val.size()+i, PyTuple_GetItem(prof.py_args,i));
        
        PyObject *pyresult = PyObject_CallObject(prof.py_profile, myvals);
        if (pyresult == NULL) {
            ERROR("can't evaluate python function");
        }
        double cppresult = PyFloat_AsDouble(pyresult);
        Py_XDECREF(pyresult);
//        Py_XDECREF(myvals);
        return cppresult;
    }
    
    

};

#endif