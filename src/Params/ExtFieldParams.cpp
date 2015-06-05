#include "ExtFieldParams.h"

#include <cmath>

using namespace std;

ExtFieldParams::ExtFieldParams(PicParams& params, InputData &ifile) :
geometry(params.geometry),
conv_fac(params.conv_fac)
{

    // -----------------
    // ExtFields properties
    // -----------------
    unsigned int numExtFields=ifile.nComponents("ExtField");
    for (unsigned int n_extfield = 0; n_extfield < numExtFields; n_extfield++) {
        ExtFieldStructure tmpExtField;
        ifile.extract("field",tmpExtField.fields,"ExtField",n_extfield);
        ifile.extract("profile",tmpExtField.profile,"ExtField",n_extfield);
        
        if (tmpExtField.profile.empty()) {
            PyObject *mypy = ifile.extract_py("profile","ExtField",n_extfield);
            if (mypy && PyCallable_Check(mypy)) {
                tmpExtField.py_profile=mypy;
                tmpExtField.profile="python";
            } else{
                ERROR(" ExtField #"<<n_extfield<<": parameter 'profile' not understood");
            }
        } else {
            if( !ifile.extract("magnitude",tmpExtField.magnitude,"ExtField",n_extfield) ) {
                ERROR("ExtField #"<<n_extfield<<": parameter 'magnitude' not provided");
            }
            ifile.extract("int_params",tmpExtField.int_params,"ExtField",n_extfield);
            ifile.extract("double_params",tmpExtField.double_params,"ExtField",n_extfield);
            ifile.extract("length_params_x",tmpExtField.length_params_x,"ExtField",n_extfield);
            ifile.extract("length_params_y",tmpExtField.length_params_y,"ExtField",n_extfield);
            ifile.extract("length_params_z",tmpExtField.length_params_z,"ExtField",n_extfield);
            
            transform(tmpExtField.length_params_x.begin(), tmpExtField.length_params_x.end(), tmpExtField.length_params_x.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
            transform(tmpExtField.length_params_y.begin(), tmpExtField.length_params_y.end(), tmpExtField.length_params_y.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
            transform(tmpExtField.length_params_z.begin(), tmpExtField.length_params_z.end(), tmpExtField.length_params_z.begin(),
                      bind1st(multiplies<double>(),params.conv_fac));
        }
        structs.push_back(tmpExtField);
    }
    
}

