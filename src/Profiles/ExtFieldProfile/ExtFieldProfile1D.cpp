#include "ExtFieldProfile1D.h"


ExtFieldProfile1D::ExtFieldProfile1D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct){

    //checking for errors
    if (my_struct.profile == "constant") {
        if (my_struct.double_params.size()<1)
            ERROR("double params size wrong " );

    } else if (my_struct.profile == "magexpansion") {
        // ---------------------------------------------------------------------------------
        // Charles magnetic field profile for Liang simulations
        // Top-hat profile :
        // int_params[0]      = order of the Gaussian
        // double_params[0]   = background Bfield amplitude
        // double_params[1]   = maximum Bfield amplitude
        // length_params_x[0] = position of the maximum of the B-field
        // length_params_x[1] = FWHM of the magnetic field (Gaussian) distribution
        // ---------------------------------------------------------------------------------
        if (my_struct.double_params.size()<1)
            ERROR("one int_params must be defined for Charles profile" );
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
    
    else if (my_struct.profile == "magexpansion") {
        // ---------------------------------------------------------------------------------
        // Charles magnetic field profile for Liang simulations
        // Top-hat profile :
        // int_params[0]      = order of the Gaussian
        // double_params[0]   = background Bfield amplitude
        // double_params[1]   = maximum Bfield amplitude
        // length_params_x[0] = position of the maximum of the B-field
        // length_params_x[1] = FWHM of the magnetic field (Gaussian) distribution
        // ---------------------------------------------------------------------------------
        int    N     = my_struct.int_params[0];
        double B0    = my_struct.double_params[0];
        double Bmax  = my_struct.double_params[1];
        double x0    = my_struct.length_params_x[0];
        double L     = my_struct.length_params_x[1];
        double sigma = pow(L/2,N)/log(2.0);
        double x     = x_cell[0]-x0;
        
        return B0 + Bmax * exp(-pow(x,N)/sigma);
            
    }
    
    else {
        return 0;
    }
};
