#include "ExtFieldProfile2D.h"

using namespace std;

ExtFieldProfile2D::ExtFieldProfile2D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct) {
    //checking for errors
    if (my_struct.profile == "constant") {
        if (my_struct.double_params.size()<1) 
            ERROR("double params size wrong " );
    } else {
        ERROR("unknown or empty profile :" << my_struct.profile );
    }
    
}


double ExtFieldProfile2D::operator() (vector<double> x_cell) {

    if (my_struct.profile == "constant") {
        return my_struct.double_params[0];
    }

    return 0;
}

