#include "InputData.h"
#include <sstream>
#include <vector>

extern "C" {
#include "pyinit.h"
#include "pycontrol.h"
}

using namespace std;

InputData::InputData(SmileiMPI *smpi, std::vector<std::string> namelistsFiles):
namelist(""),
py_namelist(NULL)
{
    Py_Initialize();
    
    py_namelist = PyImport_AddModule("__main__");
    
    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(py_namelist, "smilei_mpi_rank", smpi->getRank());
    
    // Running pyinit.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyinit_py), Python_pyinit_py_len), "pyinit.py");
    
    // Running the namelists
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
    
    // Running pycontrol.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pycontrol_py), Python_pycontrol_py_len),"pycontrol.py");
    
    // Now the string "namelist" contains all the python files concatenated
    // It is written as a file, by default "smilei.py"
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
    MESSAGE(1,"Passing through python " << name);
    DEBUG(">>>>>>>>>>>>>>> passing this to python:\n" <<command);
    int retval=PyRun_SimpleString(command.c_str());
    DEBUG("<<<<<<<<<<<<<<< from " << name);
    if (retval==-1) {
        ERROR("error parsing "<< name);
    }
}

//! retrieve python object
PyObject* InputData::extract_py(string name, string component, int nComponent) {    
//    DEBUG("[" << name << "] [" << component << "]");
    if (name.find(" ")!= string::npos || component.find(" ")!= string::npos) {
        WARNING("asking for [" << name << "] [" << component << "] : it has white inside: please fix the code");
    }
    
    PyObject *py_obj=py_namelist;
    string list = "list";
    // If component requested
    if (!component.empty()) {
        // Get the selected component (e.g. "Species" or "Laser")
        py_obj = PyObject_GetAttrString(py_namelist,component.c_str());
        PyTools::checkPyError();
        // Error if not found
        if (!py_obj) ERROR("Component "<<component<<" not found in namelist");
        // Get the "list" that contains the list of the objects in this component
        py_obj = PyObject_GetAttrString(py_obj,list.c_str());
        // If list successfully found
        if (PyList_Check(py_obj) || PyTuple_Check(py_obj)) {
            int len = PySequence_Size(py_obj);
            if (len > 0) { 
                if (len >= nComponent) {
                    PyObject* seq = PySequence_Fast(py_obj, "expected a sequence");
                    py_obj = PySequence_Fast_GET_ITEM(seq, nComponent);
                    Py_DECREF(seq);
                } else {
                    ERROR("component " << component << " is not big enough");
                }
            }
        } else {
            py_obj=NULL;
        }
    }
    PyObject *py_return=PyObject_GetAttrString(py_obj,name.c_str());
    PyTools::checkPyError();
    return py_return;

}

//! retrieve a vector of python objects
vector<PyObject*> InputData::extract_pyVec(string name, string component, int nComponent) {
    vector<PyObject*> retvec;
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

bool InputData::existComponent(std::string component, unsigned int nComponent) {
    if (component.find(" ")!= string::npos) {
        ERROR("[" << component << "] has whitespace inside: please fix the code");
    }
    // Get the selected component (e.g. "Species" or "Laser")
    string list = "list";
    PyObject *py_obj = PyObject_GetAttrString(py_namelist,component.c_str());
    PyTools::checkPyError();
    if (py_obj) {
        // Get the "list" that contains the list of the objects in this component
        py_obj = PyObject_GetAttrString(py_obj,list.c_str());
        if (py_obj) {
            if (PyList_Check(py_obj)) {
                if (PySequence_Size(py_obj) > nComponent) {
                    return true;
                }
            }
        }
    }
    return false;
}



