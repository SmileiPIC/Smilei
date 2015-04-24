#include "VelocityProfile1D.h"
#include "Tools.h"


VelocityProfile1D::VelocityProfile1D(ProfileSpecies &my_prof_params) : VelocityProfile(my_prof_params) {

    // -------------------------
    // Constant velocity profile
    // -------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density (default is sim_length-vacuum_length[0])
    if (prof_params.profile=="constant") {
        // nothing done here, by default: vacuum_length[0] = 0, dens_length_x[0] = 0
    }
    // ---------------------------------------------------------------------------------
    // Charles magnetic field profile for Liang simulations
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_x[0] = position of the maximum magnetic field
    // length_params_x[1] = Length of the magnetic gradient
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile == "magexpansion") {
        //if (prof_params.int_params.size()<2)
        //    ERROR("two int_params must be defined for Charles velocity profile" );
        if (prof_params.double_params.size()<4)
            ERROR("four double_params must be defined for Charles velocity profile" );
        if (prof_params.length_params_x.size()<2)
            ERROR("two length_params_x must be defined for Charles velocity profile" );
    } 
}

double VelocityProfile1D::operator() (std::vector<double> x_cell) {
    
    // ------------------------
    // Constant density profile
    // ------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density (default is sim_length-vacuum_length[0])
    if (prof_params.profile=="constant") {
        if (   (x_cell[0]>prof_params.vacuum_length[0])
            && (x_cell[0]<prof_params.vacuum_length[0]+prof_params.length_params_x[0]) ) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    
    // ------------------------
    // Charles velocity profile
    // ------------------------
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_x[0] = position of the maximum magnetic field
    // length_params_x[1] = Length of the magnetic gradient
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile=="magexpansion") {
        //int    N     = prof_params.int_params[0];
        //int    m     = prof_params.int_params[1];
        double n0    = prof_params.double_params[0];
        double P0    = prof_params.double_params[1];
        double B0    = prof_params.double_params[2];
        double Bmax  = prof_params.double_params[3];
        //double theta = prof_params.double_params[3];
        double x0    = prof_params.length_params_x[0];
        double L     = prof_params.length_params_x[1];
        //double alpha = Bmax/nmax * (double)(N)/sigma;
        //double sigma = pow(L/2.0,N)/log(2.0);
        double x     = x_cell[0]-x0;
	double tiny  = 1e-10*L;
	if (Bmax == 0.) {
		double Bm = sqrt(pow(B0,2) + 2*P0)-B0;
		double B  = B0 + Bm/pow(cosh(x/L),2);
		double A  = B0*x + Bm*L*tanh(x/L);
		double DP = P0 + pow(B0,2)/2 - pow(B,2)/2;
		//if (abs(x)<tiny) {     // X=0 -> velocity is 0 imposed here to avoid 0/0
		//	return (0);
		//}
		//else {	
		double v     = -2*Bm/(L*n0)*tanh(x/L) /(pow(cosh(x/L),2))*exp( 2*A*Bm/L*tanh(x/L) /(DP*pow(cosh(x/L),2)) );
        	if (abs(v)>1.0) ERROR("Velocity profile exceeding c");
		return v;
		//}
	}
	else {	
		double Bm = Bmax;
		double B  = B0 + Bm/pow(cosh(x/L),2);
		double A  = B0*x + Bm*L*tanh(x/L);
		double DP = P0 + pow(B0,2)/2 - pow(B,2)/2;
		//if (abs(x)<tiny) {
		//	return (0);    // X=0 -> velocity is 0 imposed here to avoid 0/0
		//}
		//else {	
		double v     = -2*Bm/(L*n0)*tanh(x/L) /(pow(cosh(x/L),2))*exp( 2*A*Bm/L*tanh(x/L) /(DP*pow(cosh(x/L),2)) );
        	if (abs(v)>1.0) ERROR("Velocity profile exceeding c");
		return v;
		//}
	}
    }
    
    return 1;
};
