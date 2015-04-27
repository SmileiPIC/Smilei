#include "InputData.h"
#include <sstream>
#include <vector>

extern "C" {
#include "pyinit.h"
}

using namespace std;

InputData::InputData():namelist("") {
    Py_Initialize();
    PyRun_SimpleString(reinterpret_cast<const char*>(Python_pyinit_py));
}

InputData::~InputData() {
    Py_Finalize();
}


bool BothAreSpaces(char lhs, char rhs) {
    return (lhs == rhs) && (lhs == ' ');
}

string InputData::cleanString(string str) {
    str=str.substr(0, str.find('#'));
    str=str.substr(0, str.find('!'));
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    const string whiteSpaces( " \f\n\r\t\v" );
    size_t pos = str.find_last_not_of( whiteSpaces );
    str.erase( pos + 1 );
    pos = str.find_first_not_of( whiteSpaces );
    str.erase( 0, pos );
    string::iterator new_end = unique(str.begin(), str.end(), BothAreSpaces);
    str.erase(new_end, str.end());
    return str;
}

void InputData::write(string filename=string()) {
    ofstream ostr(filename.c_str());
    if (ostr.is_open()) {
        ostr << "# smilei " << __VERSION << endl << endl;
        write(ostr);
    } else {
        write(cerr);
    }
}

void InputData::write(ostream &ostr) {
    for(vector<pair<string , vector<pair<string,string> > > >::iterator it_type = allData.begin(); it_type != allData.end(); it_type++) {
        if (!it_type->first.empty()) ostr << it_type->first << endl;
        for(vector<pair<string, string> >::iterator it_type2 = it_type->second.begin(); it_type2 != it_type->second.end(); it_type2++) {
            if (!it_type->first.empty()) ostr << "\t";
            ostr << it_type2->first << " = " << it_type2->second << endl;
        }
        if (!it_type->first.empty()) ostr << "end" << endl;
        ostr << endl;
    }
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
    PyObject* py_val = PyObject_GetAttrString(py_namelist,name.c_str());
    return py_val;
}    

//! retrieve python object
vector<PyObject*> InputData::py_vec_from_string(string name, string group, int occurrenceItem, int occurrenceGroup) {
    vector<PyObject*> retvec;
    PyObject* py_val = PyObject_GetAttrString(py_namelist,name.c_str());
    
    if (py_val) {        
        if (!PyTuple_Check(py_val)) {
            retvec.push_back(py_val);
        } else {
            PyObject* seq = PySequence_Fast(py_val, "expected a sequence");
            int len = PySequence_Size(py_val);
            retvec.resize(len);
            for (int i = 0; i < len; i++) {
                PyObject* item = PySequence_Fast_GET_ITEM(seq, i);
                retvec[i]=item;
                DEBUG(" here " << i);
            }
            Py_DECREF(seq);
        }      
    }    
    return retvec;
}    


void InputData::parseStream() {
    
    //! we let python execute the namelist
    int retval=PyRun_SimpleString(namelist.c_str());
    if (retval==-1) {
        ERROR("error parsing namelist")
    }
    
    //this is  apython function described in pyinit.py
    PyRun_SimpleString("check_namelist()");
    
    // we store in a pyobject the smilei class of the namelist
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

bool InputData::existGroup(string groupName, unsigned int occurrenceGroup) {
    unsigned int n_occur_group=0;
    for (vector<pair<string , vector<pair<string,string> > > >::iterator  it_type = allData.begin(); it_type != allData.end(); it_type++) {
        if (groupName == it_type->first) {
            if (occurrenceGroup==n_occur_group) {
                return true;
            }
            n_occur_group++;
        }
    }
    return false;
}


