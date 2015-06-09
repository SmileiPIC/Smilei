#include "InputData.h"
#include <sstream>
#include <vector>

#include "pyinit.h"
#include "pyprofiles.h"
#include "pycontrol.h"

#include <sys/time.h>
#include <sys/resource.h>

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

    // Running pyfunctons.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyinit_py), Python_pyprofiles_py_len), "pyprofiles.py");

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

    PyTools::runPyFunction("smilei_check");
    
    
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
    if (Py_IsInitialized())
        Py_Finalize();
}

//! run script
void InputData::pyRunScript(string command, string name) {
    PyTools::checkPyError();
    namelist+=command;
    MESSAGE(1,"Passing to python " << name);
    DEBUG(">>>>>>>>>>>>>>> passing this to python:\n" <<command);
    int retval=PyRun_SimpleString(command.c_str());
    DEBUG("<<<<<<<<<<<<<<< from " << name);
    if (retval==-1) {
        ERROR("error parsing "<< name);
        PyTools::checkPyError();
    }
}

//! retrieve python object
PyObject* InputData::extract_py(string name, string component, int nComponent) {
//    DEBUG("[" << name << "] [" << component << "]");
    if (name.find(" ")!= string::npos || component.find(" ")!= string::npos) {
        WARNING("asking for [" << name << "] [" << component << "] : it has whitespace inside: please fix the code");
    }
    PyObject *py_obj=py_namelist;
    // If component requested
    if (!component.empty()) {
        // Get the selected component (e.g. "Species" or "Laser")
        py_obj = PyObject_GetAttrString(py_namelist,component.c_str());
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

int InputData::nComponents(std::string componentName) {
    // Get the selected component (e.g. "Species" or "Laser")
    PyObject *py_obj = PyObject_GetAttrString(py_namelist,componentName.c_str());
    PyTools::checkPyError();
    int retval = PyObject_Length(py_obj);
    HEREIAM(componentName << " " << retval);
    return retval;
}


//! run the python functions cleanup (user defined) and keep_python_running (in pycontrol.py)
void InputData::cleanup() {
    
    // call cleanup function from the user namelist (it can be used to free some memory 
    // from the python side) while keeping the interpreter running
    MESSAGE(1,"Checking for cleanup() function:");
    PyTools::runPyFunction("cleanup");
    // this will reset error in python in case cleanup doesn't exists
    PyErr_Clear();

    // this function is defined in the Python/pyontrol.py file and should return false if we can close
    // the python interpreter
    MESSAGE(1,"Calling python keep_python_running() :");    
    if (PyTools::runPyFunction<bool>("keep_python_running")) {
        MESSAGE(2,"Keeping Python interpreter alive");
    } else {
        MESSAGE(2,"Closing Python");
        PyErr_Print();
        Py_Finalize();
    }
}

