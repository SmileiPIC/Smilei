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
    while (ifile.existGroup(groupName,n_extfield)) {
        
        ExtFieldStructure tmpExtField;
        ifile.extract("field",tmpExtField.fields,groupName,0,n_extfield);
        ifile.extract("profile",tmpExtField.profile,groupName,0,n_extfield);
        ifile.extract("int_params",tmpExtField.int_params,groupName,0,n_extfield);
        ifile.extract("double_params",tmpExtField.double_params,groupName,0,n_extfield);
        ifile.extract("length_params_x",tmpExtField.length_params_x,groupName,0,n_extfield);
        ifile.extract("length_params_y",tmpExtField.length_params_y,groupName,0,n_extfield);
        ifile.extract("length_params_z",tmpExtField.length_params_z,groupName,0,n_extfield);

        transform(tmpExtField.length_params_x.begin(), tmpExtField.length_params_x.end(), tmpExtField.length_params_x.begin(),
                  bind1st(multiplies<double>(),params.conv_fac));
        transform(tmpExtField.length_params_y.begin(), tmpExtField.length_params_y.end(), tmpExtField.length_params_y.begin(),
                  bind1st(multiplies<double>(),params.conv_fac));
        transform(tmpExtField.length_params_z.begin(), tmpExtField.length_params_z.end(), tmpExtField.length_params_z.begin(),
                  bind1st(multiplies<double>(),params.conv_fac));
        
        structs.push_back(tmpExtField);

        n_extfield++;
    }
    
}

