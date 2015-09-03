#include "ExtFieldParams.h"

#include <cmath>

using namespace std;

ExtFieldParams::ExtFieldParams(Params& params) :
geometry(params.geometry)
{

    // -----------------
    // ExtFields properties
    // -----------------
    unsigned int numExtFields=PyTools::nComponents("ExtField");
    for (unsigned int n_extfield = 0; n_extfield < numExtFields; n_extfield++) {
        ExtFieldStructure tmpExtField;
        if( !PyTools::extract("field",tmpExtField.fields,"ExtField",n_extfield)) {
            ERROR("ExtField #"<<n_extfield<<": parameter 'field' not provided'");
        }
        
        // If profile is a float
        if( PyTools::extract("profile", tmpExtField.profile, "ExtField", n_extfield) ) {
            string xyz = "x";
            if(geometry=="2d3v") xyz = "x,y";
            if(geometry=="3d3v") xyz = "x,y,z";
            // redefine the profile as a constant function instead of float
            PyTools::checkPyError();
            ostringstream command;
            command.str("");
            command << "ExtField["<<n_extfield<<"].profile=lambda "<<xyz<<":" << tmpExtField.profile;
            if( !PyRun_SimpleString(command.str().c_str()) ) PyTools::checkPyError();
        }
        // Now import the profile as a python function
        PyObject *mypy = PyTools::extract_py("profile","ExtField",n_extfield);
        if (mypy && PyCallable_Check(mypy)) {
            tmpExtField.py_profile=mypy;
        } else{
            ERROR(" ExtField #"<<n_extfield<<": parameter 'profile' not understood");
        }
        structs.push_back(tmpExtField);
    }

    
    unsigned int numAntenna=PyTools::nComponents("Antenna");
    for (unsigned int n_antenna = 0; n_antenna < numAntenna; n_antenna++) {
        HEREIAM(n_antenna);
        
        AntennaStructure tmpProf;
        tmpProf.my_field = NULL;
        if( !PyTools::extract("field",tmpProf.field,"Antenna",n_antenna)) {
            ERROR("ExtField #"<<n_antenna<<": parameter 'field' not provided'");
        }
        if (tmpProf.field != "Jx" && tmpProf.field != "Jy" && tmpProf.field != "Jz")
            ERROR("Antenna field must be one of J{x,y,z}");
        
        // Now import the space profile as a python function
        PyObject *mypy = PyTools::extract_py("space_profile","Antenna",n_antenna);
        if (mypy && PyCallable_Check(mypy)) {
            tmpProf.space_profile.py_profile=mypy;
        } else{
            ERROR(" Antenna #"<<n_antenna<<": parameter 'space_profile' not understood");
        }

        // Now import the time profile as a python function
        PyObject *mytpy = PyTools::extract_py("time_profile","Antenna",n_antenna);
        if (mytpy && PyCallable_Check(mytpy)) {
            tmpProf.time_profile.py_profile=mytpy;
        } else{
            ERROR(" Antenna #"<<n_antenna<<": parameter 'time_profile' not understood");
        }
        
        antennas.push_back(tmpProf);
    }

}

