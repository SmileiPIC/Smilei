#include "InputData.h"
#include <sstream>
#include <vector>

extern "C" {
#include "pyinit.h"
}

using namespace std;

InputData::InputData():namelist(""), py_namelist(NULL) {
    Py_Initialize();
    // this will read the file pyinit.h
    PyRun_SimpleString(reinterpret_cast<const char*>(Python_pyinit_py));
}

InputData::~InputData() {
    Py_Finalize();
}

//! get bool from python
bool InputData::extract(string name, bool &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    DEBUG(name << " " << py_val);
    if(py_val) {
        if (PyBool_Check(py_val)) {
            return py_val==Py_True;
        } else {
            ERROR(name << " is not a boolean");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, short int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (short int) lval;
            DEBUG(name << " : " << val);
            return true;
        } else {
            ERROR(name << " is not a short int");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, unsigned int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (unsigned int) lval;
            DEBUG(name << " : " << val);
            return true;
        } else {
            ERROR(name << " is not a unsigned int");
        }
    }
    return false;
}

//! get int from python
bool InputData::extract(string name, int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (int) lval;
            DEBUG(name << " : " << val);
            return true;
        } else {
            ERROR(name << " is not a int");
        }
    }
    return false;
}

//! get double from python
bool InputData::extract(string name, double &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyFloat_Check(py_val)) {
            val = PyFloat_AsDouble(py_val);
            DEBUG(name << " : " << val);
            return true;
        } else {
            ERROR(name << " is not a double");
        }
    }
    return false;
}

//! get string from python
bool InputData::extract(string name, string &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyString_Check(py_val)) {
            const char* s = PyString_AsString(py_val);
            val=string(s);
            DEBUG(name << " : " << val);
            return true;
        } else {
            ERROR(name << " is not a string");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, vector<unsigned int> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=py_vec_from_string(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyInt_Check(pyvec[i])) {
            long int lval = PyInt_AsLong(pyvec[i]);
            val[i] = (unsigned int) lval;
        } else {
            ERROR("reading float in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get int from python
bool InputData::extract(string name, vector<int> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=py_vec_from_string(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyInt_Check(pyvec[i])) {
            long int lval = PyInt_AsLong(pyvec[i]);
            val[i] = (int) lval;
        } else {
            ERROR("reading float in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get double from python
bool InputData::extract(string name, vector<double> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=py_vec_from_string(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyFloat_Check(pyvec[i])) {
            val[i] = PyFloat_AsDouble(pyvec[i]);
        } else {
            ERROR("reading float in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get string from python
bool InputData::extract(string name, vector<string> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=py_vec_from_string(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyString_Check(pyvec[i])) {
            val[i] = PyFloat_AsDouble(pyvec[i]);
        } else {
            ERROR("reading string in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! retrieve python object
PyObject* InputData::py_val_from_string(string name, string group, int occurrenceItem, int occurrenceGroup) {    
    PyObject *py_obj=py_namelist;
    if (!group.empty()) {
        py_obj = PyObject_GetAttrString(py_namelist,group.c_str());
        if (py_obj) {
            if (PyList_Check(py_obj)) {
                int len = PySequence_Size(py_obj);
                if (len >= occurrenceGroup) {
                    PyObject* seq = PySequence_Fast(py_obj, "expected a sequence");
                    py_obj = PySequence_Fast_GET_ITEM(seq, occurrenceGroup);
                    Py_DECREF(seq);
                } else {
                    ERROR("group " << group << " is not big enough");
                }
            }
        } else {
            ERROR("group " << group << " does not exists");
        }
    }
    PyObject* py_val=NULL;
    if (py_obj) {
        py_val = PyObject_GetAttrString(py_obj,name.c_str());
        if (occurrenceItem>0) {
            if (PyList_Check(py_val)) {
                int len = PySequence_Size(py_val);
                if (len >= occurrenceGroup) {
                    PyObject* seq = PySequence_Fast(py_val, "expected a sequence");
                    py_val = PySequence_Fast_GET_ITEM(seq, occurrenceItem);
                    Py_DECREF(seq);
                } else {
                    ERROR(name << " is not big enough");
                }
            } else {
                ERROR(name << " is not a list");
            }

        }
    } else {
        ERROR("something bad happened " << name << "in group " << group);
    }
    return py_val;
}    

//! retrieve a vector of python objects
vector<PyObject*> InputData::py_vec_from_string(string name, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = py_val_from_string(name,group,occurrenceItem,occurrenceGroup);
    vector<PyObject*> retvec;
    if (py_val) {      
        if (!PyTuple_Check(py_val)) {
            retvec.push_back(py_val);
            WARNING(name << " should be a tuple, not a scalar : fix it");
        } else {
            PyObject* seq = PySequence_Fast(py_val, "expected a sequence");
            int len = PySequence_Size(py_val);
            retvec.resize(len);
            for (int i = 0; i < len; i++) {
                PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
                retvec[i]=item;
            }
            Py_DECREF(seq);
        }      
    }    
    return retvec;
}    


bool InputData::existGroup(std::string group, unsigned int occurrenceGroup) {
    PyObject *py_obj = PyObject_GetAttrString(py_namelist,group.c_str());
    if (py_obj) {
        if (PyList_Check(py_obj)) {
            int len = PySequence_Size(py_obj);
            if (len > occurrenceGroup) {
                DEBUG("here " << len << " " << occurrenceGroup);
                return true;
            }
        }
    }
    return false;
}

void InputData::parseStream() {
    
    //! we let python execute the namelist
    int retval=PyRun_SimpleString(namelist.c_str());
    if (retval==-1) {
        ERROR("error parsing namelist")
    }
    
    PyObject* myFunction = PyObject_GetAttrString(PyImport_AddModule("__main__"),(char*)"get_smilei");
    py_namelist = PyObject_CallFunction(myFunction,const_cast<char *>(""));
    
    if (!py_namelist) {
        ERROR("no smilei class defined")
    }
}


void InputData::readFile(string filename) {
    
    ifstream istr(filename.c_str());
    
    string strLine ="";
    namelist.clear();
    
    if (istr.is_open()) {
        while (getline(istr, strLine)) {
            namelist += strLine + "\n";
        }
    } else {
        ERROR("File " << filename << " does not exists");
    }
    namelist +="\n";    
}



