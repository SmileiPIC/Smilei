#include "InputData.h"
#include <sstream>
#include <vector>

#include "pyinit.pyh"
#include "pyprofiles.pyh"
#include "pycontrol.pyh"

#include <sys/time.h>
#include <sys/resource.h>

using namespace std;

InputData::InputData(SmileiMPI *smpi, std::vector<std::string> namelistsFiles):
namelist("")
{
    if (!Py_IsInitialized()) Py_Initialize();
    
    // here we add the rank, in case some script need it
    PyModule_AddIntConstant(PyImport_AddModule("__main__"), "smilei_mpi_rank", smpi->getRank());
    
    // First, we tell python to filter the ctrl-C kill command (or it would prevent to kill the code execution).
    // This is done separately from other scripts because we don't want it in the concatenated python namelist.
    PyTools::checkPyError();
    string command = "import signal\nsignal.signal(signal.SIGINT, signal.SIG_DFL)";
    if( !PyRun_SimpleString(command.c_str()) ) PyTools::checkPyError();
    
    // Running pyinit.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyinit_py), Python_pyinit_py_len), "pyinit.py");
    
    // Running pyfunctons.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pyprofiles_py), Python_pyprofiles_py_len), "pyprofiles.py");
    
    // Running the namelists
    pyRunScript("############### BEGIN USER NAMELISTS ###############\n");
    for (vector<string>::iterator it=namelistsFiles.begin(); it!=namelistsFiles.end(); it++) {
        MESSAGE(1,"Reading file " << *it);
        string strNamelist="";
        if (smpi->isMaster()) {
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
    pyRunScript("################ END USER NAMELISTS ################\n");    
    // Running pycontrol.py
    pyRunScript(string(reinterpret_cast<const char*>(Python_pycontrol_py), Python_pycontrol_py_len),"pycontrol.py");
    
    PyTools::runPyFunction("_smilei_check");
    
    
    // Now the string "namelist" contains all the python files concatenated
    // It is written as a file: smilei.py
    if (smpi->isMaster()) {
        ofstream out_namelist("smilei.py");
        if (out_namelist.is_open()) {
            out_namelist << namelist;
            out_namelist.close();
        }
    }
}

InputData::~InputData() {
    if (Py_IsInitialized())
        Py_Finalize();
}

//! Run string as python script and add to namelist
void InputData::pyRunScript(string command, string name) {
    PyTools::checkPyError();
    namelist+=command;
    if (name.size()>0)  MESSAGE(1,"Passing to python " << name);
    DEBUG(">>>>>>>>>>>>>>> passing this to python:\n" <<command);
    int retval=PyRun_SimpleString(command.c_str());
    DEBUG("<<<<<<<<<<<<<<< from " << name);
    if (retval==-1) {
        ERROR("error parsing "<< name);
        PyTools::checkPyError();
    }
}

//! run the python functions cleanup (user defined) and _keep_python_running (in pycontrol.py)
void InputData::cleanup() {
    
    // call cleanup function from the user namelist (it can be used to free some memory 
    // from the python side) while keeping the interpreter running
    MESSAGE(1,"Checking for cleanup() function:");
    PyTools::runPyFunction("cleanup");
    // this will reset error in python in case cleanup doesn't exists
    PyErr_Clear();
    
    // this function is defined in the Python/pyontrol.py file and should return false if we can close
    // the python interpreter
    MESSAGE(1,"Calling python _keep_python_running() :");    
    if (PyTools::runPyFunction<bool>("_keep_python_running")) {
        MESSAGE(2,"Keeping Python interpreter alive");
    } else {
        MESSAGE(2,"Closing Python");
        PyErr_Print();
        Py_Finalize();
    }
}

