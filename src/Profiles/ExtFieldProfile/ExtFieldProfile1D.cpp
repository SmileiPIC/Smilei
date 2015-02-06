#include "ExtFieldProfile1D.h"


ExtFieldProfile1D::ExtFieldProfile1D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct){

    //checking for errors
    if (my_struct.profile == "constant") {
        if (my_struct.double_params.size()<1)
            ERROR("double params size wrong " );

    } else if (my_struct.profile == "charles") {
        // Top-Hat profile (Charles Ruyer simulations)
        if (my_struct.double_params.size()<2)
            ERROR("two double_params must be defined for Charles profile" );
        if (my_struct.length_params_x.size()<2)
            ERROR("three length_params_x must be defined for Charles profile" );
        

    } else {
        ERROR("unknown or empty profile :" << my_struct.profile );
    }


}


double ExtFieldProfile1D::operator() (std::vector<double> x_cell) {

    if (my_struct.profile == "constant") {
        return my_struct.double_params[0];
    }
    
    else if (my_struct.profile == "charles") {
        // ---------------------------------------------------------------------------------
        // Charles magnetic field profile for Liang simulations
        // Top-hat profile :
        // double_params[0]   = background Bfield amplitude
        // double_params[1]   = maximum Bfield amplitude
        // length_params_x[0] = length of the first plateau
        // length_params_x[1] = length of the 2nd plateau (where the field is cst & maximum)
        // ---------------------------------------------------------------------------------
        double B0   = my_struct.double_params[0];
        double Bmax = my_struct.double_params[1];
        double x0 = my_struct.length_params_x[0];
        double x1 = my_struct.length_params_x[1];
        
        if ( (x0<x_cell[0]) && (x_cell[0]<x1) ) {
            return Bmax;
        } else {
            return B0;
        }
        
    } else {
        return 0;
    }
};
