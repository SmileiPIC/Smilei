#include "InputData.h"
#include <sstream>
#include <vector>

extern "C" {
#include "pyinit.h"
#include "pycontrol.h"
}

using namespace std;

InputData::InputData(SmileiMPI *smpi, std::vector<std::string> namelistsFiles): namelist(""), py_namelist(NULL) {
    Py_Initialize();

    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_rank", smpi->getRank());
    
    // here we add the main python class definitions from pyinit.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyinit_py), Python_pyinit_py_len), "pyinit.py");
    
    PyObject *rank=Py_BuildValue("i", smpi->getRank());
    
    
    // we read the file(s)
    for (vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++) {
        MESSAGE("Reading file " << *it);
        string strNamelist="";
        if (smpi->isMaster()) {
            HEREIAM("");
            ifstream istr(it->c_str());        
            if (istr.is_open()) {
                string oneLine;
                while (getline(istr, oneLine)) {
                    strNamelist += oneLine + "\n";
                }
            } else {
                ERROR("File " << (*it) << " does not exists");
            }
            strNamelist +="\n";
        }
        smpi->bcast(strNamelist);
        pyRunScript(strNamelist,(*it));
    }
    
    
    // here we add the check namelist stuff from pycontrol.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pycontrol_py), Python_pycontrol_py_len),"pycontrol.py");

    
    

    int retval=PyRun_SimpleString(namelist.c_str());
    if (retval==-1) {
        ERROR("error parsing namelist")
    }

    PyObject* myFunction = PyObject_GetAttrString(PyImport_AddModule("__main__"),(char*)"get_smilei");
    py_namelist = PyObject_CallFunction(myFunction,const_cast<char *>(""));
    if (!py_namelist) {
        ERROR("no smilei class defined, but we should never get here...");
    }
    
    if (smpi->isMaster()) {
        string file_namelist_out="smilei.py";
        extract("output_script", file_namelist_out);
        
        ofstream out(file_namelist_out.c_str());
        out << namelist;
        out.close();
    }        
}

InputData::~InputData() {
    Py_Finalize();
}

//! run script
void InputData::pyRunScript(string command, string name) {
    namelist+=command;
    MESSAGE(1,"passing to python " << name);
    DEBUG(">>>>>>>>>>>>>>> passing this to python:\n" <<command);
    int retval=PyRun_SimpleString(command.c_str());
    DEBUG("<<<<<<<<<<<<<<< from " << name);
    if (retval==-1) {
        ERROR("error parsing "<< name);
    }
}

//! get bool from python
bool InputData::extract(string name, bool &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if(py_val) {
        if (PyBool_Check(py_val)) {
            val=(py_val==Py_True);
            HEREIAM(name << " " << val);
            return true;
        } else {
            DEBUG(name << " is not a boolean");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, short int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (short int) lval;
            return true;
        } else {
            DEBUG(name << " is not a short int");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, unsigned int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (unsigned int) lval;
            return true;
        } else {
            DEBUG(name << " is not a unsigned int");
        }
    }
    return false;
}

//! get int from python
bool InputData::extract(string name, int &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyInt_Check(py_val)) {
            long int lval = PyInt_AsLong(py_val);
            val = (int) lval;
            return true;
        } else {
            DEBUG(name << " is not a int");
        }
    }
    return false;
}

//! get double from python
bool InputData::extract(string name, double &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyFloat_Check(py_val)) {
            val = PyFloat_AsDouble(py_val);
            return true;
        } else if (PyInt_Check(py_val)) {
            val=(double) PyInt_AsLong(py_val);
            return true;
        } else {
            DEBUG(name << " is not a double");
        }
    }
    return false;
}

//! get string from python
bool InputData::extract(string name, string &val, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
    if (py_val) {
        if (PyString_Check(py_val)) {
            const char* s = PyString_AsString(py_val);
            val=string(s);
            return true;
        } else {
            DEBUG(name << " is not a string");
        }
    }
    return false;
}

//! get uint from python
bool InputData::extract(string name, vector<unsigned int> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=extract_pyVvec(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyInt_Check(pyvec[i])) {
            long int lval = PyInt_AsLong(pyvec[i]);
            val[i] = (unsigned int) lval;
        } else {
            DEBUG("reading unsigned int in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get int from python
bool InputData::extract(string name, vector<int> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=extract_pyVvec(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyInt_Check(pyvec[i])) {
            long int lval = PyInt_AsLong(pyvec[i]);
            val[i] = (int) lval;
        } else {
            DEBUG("reading int in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get double from python
bool InputData::extract(string name, vector<double> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=extract_pyVvec(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyFloat_Check(pyvec[i])) {
            val[i] = PyFloat_AsDouble(pyvec[i]);
        } else if (PyInt_Check(pyvec[i])) {
            val[i] = (double) PyInt_AsLong(pyvec[i]);
        } else {
            DEBUG("reading float in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! get string from python
bool InputData::extract(string name, vector<string> &val, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> pyvec=extract_pyVvec(name,group,occurrenceItem,occurrenceGroup);
    val.resize(pyvec.size());
    for (unsigned int i=0;i<pyvec.size();i++) {
        if (PyString_Check(pyvec[i])) {
            val[i]=string(PyString_AsString(pyvec[i]));
        } else {
            DEBUG("reading string in " << name << " at pos " <<i );
        }
    }
    return false;
}

//! retrieve python object
PyObject* InputData::extract_py(string name, string group, int occurrenceItem, int occurrenceGroup) {    
//    DEBUG("[" << name << "] [" << group << "]");
    if (name.find(" ")!= string::npos || group.find(" ")!= string::npos) {
        WARNING("asking for [" << name << "] [" << group << "] : it has white inside: please fix the code");
    }

    PyObject *py_obj=py_namelist;
    if (!group.empty()) {
        py_obj = PyObject_GetAttrString(py_namelist,group.c_str());
        if (py_obj) {
            if (PyList_Check(py_obj)) {
                int len = PySequence_Size(py_obj);
                if (len > 0) { 
                    if (len >= occurrenceGroup) {
                        PyObject* seq = PySequence_Fast(py_obj, "expected a sequence");
                        py_obj = PySequence_Fast_GET_ITEM(seq, occurrenceGroup);
                        Py_DECREF(seq);
                    } else {
                        ERROR("group " << group << " is not big enough");
                    }
                }
            } else {
                py_obj=NULL;
            }

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
    }
    return py_val;
}    

//! retrieve a vector of python objects
vector<PyObject*> InputData::extract_pyVvec(string name, string group, int occurrenceItem, int occurrenceGroup) {
    PyObject* py_val = extract_py(name,group,occurrenceItem,occurrenceGroup);
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
    if (group.find(" ")!= string::npos) {
        ERROR("[" << group << "] has white inside: please fix the code");
    }
    PyObject *py_obj = PyObject_GetAttrString(py_namelist,group.c_str());
    if (py_obj) {
        if (PyList_Check(py_obj)) {
            if (PySequence_Size(py_obj) > occurrenceGroup) {
                return true;
            }
        }
    }
    return false;
}



