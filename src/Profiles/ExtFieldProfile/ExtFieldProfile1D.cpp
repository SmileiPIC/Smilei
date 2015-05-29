#include "ExtFieldProfile1D.h"

using namespace std;

ExtFieldProfile1D::ExtFieldProfile1D(ExtFieldStructure &extfield_struct) : ExtFieldProfile(extfield_struct){

    //checking for errors
    if (my_struct.profile == "constant") {
        if (my_struct.double_params.size()<1)
            ERROR("double params size wrong " );
    } 
    else if (my_struct.profile == "magexpansion") {
        // ---------------------------------------------------------------------------------
        // Charles magnetic field profile for Liang simulations
        // double_params[0]   = background Bfield amplitude
        // double_params[1]   = maximum Bfield amplitude
        // double_params[2]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
        // length_params_x[0] = position of the maximum of the B-field
        // length_params_x[1] = Length of the magnetic gradient
        // ---------------------------------------------------------------------------------
        //if (my_struct.int_params.size()<1)
        //    ERROR("one int_params must be defined for Charles profile" );
        if (my_struct.double_params.size()<3)
            ERROR("three double_params must be defined for Charles profile" );
        if (my_struct.length_params_x.size()<2)
            ERROR("three length_params_x must be defined for Charles profile" );

    } 
    else if (my_struct.profile=="python") {
    }
    else {
        ERROR("unknown or empty profile :" << my_struct.profile );
    }


}


double ExtFieldProfile1D::operator() (vector<double> x_cell) {

    if (my_struct.profile == "constant") {
        return my_struct.double_params[0];
    }
    
    else if (my_struct.profile == "magexpansion") {
        // ---------------------------------------------------------------------------------
        // Charles magnetic field profile for Liang simulations
        // double_params[0]   = background Bfield amplitude
        // double_params[1]   = maximum Bfield amplitude
        // double_params[2]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
        // length_params_x[0] = position of the maximum of the B-field
        // length_params_x[1] = Length of the magnetic gradient
        // ---------------------------------------------------------------------------------
        //int    N     = my_struct.int_params[0];
        double B0    = my_struct.double_params[0];
        double Bmax  = my_struct.double_params[1];
	double P0    = my_struct.double_params[2];
        double x0    = my_struct.length_params_x[0];
        double L     = my_struct.length_params_x[1];
        double x     = x_cell[0]-x0;
	//double Bm    = Bmax;
        //double sigma = pow(L/2.0,N)/log(2.0);
	if (Bmax == 0.) {
		//double Bm = sqrt(B0*B0+2*P0)-B0;
		return B0 + (sqrt(B0*B0+2*P0)-B0)/pow(cosh(x/L),2);
	}else {	return B0 + Bmax/pow(cosh(x/L),2);
	}
        
        //return B0 + Bmax * exp(-pow(x,N)/sigma);
            
    }
    else if (my_struct.profile=="python") {
        return PyTools::runPyFunction(my_struct.py_profile, x_cell[0]);
    }
    else {
        return 0;
    }
};
