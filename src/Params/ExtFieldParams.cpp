#include "ExtFieldParams.h"

#include <cmath>

#include "Tools.h"

using namespace std;

ExtFieldParams::ExtFieldParams(PicParams& params, InputData &ifile, string groupName) :
ProfileParams(params)
{
	
    // -----------------
    // ExtFields properties
    // -----------------
    int n_extfield=0;
    while (ifile.existComponent(groupName,n_extfield)) {
        
        ExtFieldStructure tmpExtField;
        ifile.extract("field",tmpExtField.fields,groupName,n_extfield);
        ifile.extract("profile",tmpExtField.profile,groupName,n_extfield);
        
        if (tmpExtField.profile.empty()) {
            PyObject *mypy = ifile.extract_py("profile",groupName,n_extfield);
            if (mypy && PyCallable_Check(mypy)) {
                tmpExtField.py_profile=mypy;
                tmpExtField.profile="python";
            }
        } else {
            ifile.extract("int_params",tmpExtField.int_params,groupName,n_extfield);
            ifile.extract("double_params",tmpExtField.double_params,groupName,n_extfield);
            ifile.extract("length_params_x",tmpExtField.length_params_x,groupName,n_extfield);
            ifile.extract("length_params_y",tmpExtField.length_params_y,groupName,n_extfield);
            ifile.extract("length_params_z",tmpExtField.length_params_z,groupName,n_extfield);
            
            transform(tmpExtField.length_params_x.begin(), tmpExtField.length_params_x.end(), tmpExtField.length_params_x.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
            transform(tmpExtField.length_params_y.begin(), tmpExtField.length_params_y.end(), tmpExtField.length_params_y.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
            transform(tmpExtField.length_params_z.begin(), tmpExtField.length_params_z.end(), tmpExtField.length_params_z.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
        }   
        structs.push_back(tmpExtField);

        n_extfield++;
    }
    
}

