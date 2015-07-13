#include "ExtFieldParams.h"

#include <cmath>

using namespace std;

ExtFieldParams::ExtFieldParams(PicParams& params, InputData &ifile) :
geometry(params.geometry)
{

    // -----------------
    // ExtFields properties
    // -----------------
    unsigned int numExtFields=ifile.nComponents("ExtField");
    for (unsigned int n_extfield = 0; n_extfield < numExtFields; n_extfield++) {
        ExtFieldStructure tmpExtField;
        if( !ifile.extract("field",tmpExtField.fields,"ExtField",n_extfield)) {
            ERROR("ExtField #"<<n_extfield<<": parameter 'field' not provided'");
        }
        
        // If profile is a float
        if( ifile.extract("profile", tmpExtField.profile, "ExtField", n_extfield) ) {
            string xyz = "x";
            if(geometry=="2d3v") xyz = "x,y";
            // redefine the profile as a constant function instead of float
            PyTools::checkPyError();
            ostringstream command;
            command.str("");
            command << "ExtField["<<n_extfield<<"].profile=lambda "<<xyz<<":" << tmpExtField.profile;
            if( !PyRun_SimpleString(command.str().c_str()) ) PyTools::checkPyError();
        }
        // Now import the profile as a python function
        PyObject *mypy = ifile.extract_py("profile","ExtField",n_extfield);
        if (mypy && PyCallable_Check(mypy)) {
            tmpExtField.py_profile=mypy;
        } else{
            ERROR(" ExtField #"<<n_extfield<<": parameter 'profile' not understood");
        }
        
        structs.push_back(tmpExtField);
    }

}

