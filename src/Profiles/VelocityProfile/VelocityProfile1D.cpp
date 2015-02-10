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
    // Top-hat profile :
    // int_params[0]      = order of the Gaussian
    // int_params[1]      = 0 for ions, 1 for electrons
    // double_params[0]   = background density
    // double_params[1]   = maximum density
    // double_params[2]   = maximum Bfield amplitude
    // double_params[3]   = electron over ion temperature ratio
    // length_params_x[0] = position of the maximum of the B-field
    // length_params_x[1] = FWHM of the magnetic field (Gaussian) distribution
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile == "magexpansion") {
        
        if (prof_params.int_params.size()<2)
            ERROR("two int_params must be defined for Charles velocity profile" );
        if (prof_params.double_params.size()<4)
            ERROR("four double_params must be defined for Charles velocity profile" );
        if (prof_params.length_params_x.size()<2)
            ERROR("three length_params_x must be defined for Charles velocity profile" );
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
    // int_params[0]      = order of the Gaussian
    // double_params[0]   = background density
    // double_params[1]   = maximum density
    // double_params[2]   = maximum magnetic field
    // double_params[3]   = ratio of electron over ion temperature
    // length_params_x[0] = position of the maximum density
    // length_params_x[1] = FWHM of the density (Gaussian) distribution
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile=="magexpansion") {
        int    N     = prof_params.int_params[0];
        int    m     = prof_params.int_params[1];
        double n0    = prof_params.double_params[0];
        double nmax  = prof_params.double_params[1];
        double Bmax  = prof_params.double_params[2];
        double theta = prof_params.double_params[3];
        double x0    = prof_params.length_params_x[0];
        double L     = prof_params.length_params_x[1];
        double sigma = pow(L/2,N)/log(2.0);
        double x     = x_cell[0]-x0;
        double alpha = pow(-theta,m)/(1.0+theta) * Bmax/nmax * (double)(N)/sigma;
        double v     = alpha * pow(x,N-1) * exp(-pow(x,N)/sigma) / (exp(-pow(x,N)/sigma)+n0/nmax);
        
        if (abs(v)>1.0) ERROR("Velocity profile exceeding c");
        return v;
    }
    
    return 1;
};
