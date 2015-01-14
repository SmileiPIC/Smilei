#include "ExtFieldParams.h"

#include <cmath>

#include "Tools.h"

using namespace std;

ExtFieldParams::ExtFieldParams(PicParams& params, InputData &ifile) :
geometry(params.geometry)
{
	
    // -----------------
    // ExtFields properties
    // -----------------
    int n_extfield=0;
    while (ifile.existGroup("extfield",n_extfield)) {
        
        ExtFieldStructure tmpExtField;
        ifile.extract("field",tmpExtField.fields,"extfield",0,n_extfield);
        ifile.extract("profile",tmpExtField.profile,"extfield",0,n_extfield);
        ifile.extract("int_params",tmpExtField.int_params,"extfield",0,n_extfield);
        ifile.extract("double_params",tmpExtField.double_params,"extfield",0,n_extfield);
        ifile.extract("length_params",tmpExtField.length_params,"extfield",0,n_extfield);

        transform(tmpExtField.length_params.begin(), tmpExtField.length_params.end(), tmpExtField.length_params.begin(),
                  bind1st(multiplies<double>(),params.conv_fac));
        
        structs.push_back(tmpExtField);

        n_extfield++;
    }
    
}

