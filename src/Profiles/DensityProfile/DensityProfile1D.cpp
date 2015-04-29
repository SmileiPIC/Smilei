#include "DensityProfile1D.h"
#include "Tools.h"


DensityProfile1D::DensityProfile1D(SpeciesStructure &params) : DensityProfile(params){

    // ------------------------
    // Constant density profile
    // ------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density (default is sim_length-vacuum_length[0])
    if (species_param.dens_profile.profile=="constant") {
        // nothing done here, by default: vacuum_length[0] = 0, dens_length_x[0] = 0
    }
    
    // ---------------------------
    // Trapezoidal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_x[1]   : length of the left slope (default value is zero)
    // dens_length_x[2]   : length of the right slope (default value is the rising slope value species_length[1])
    else if (species_param.dens_profile.profile=="trapezoidal") {
        
        if (species_param.dens_profile.length_params_x.size()<2) {
            species_param.dens_profile.length_params_x.resize(2);
            species_param.dens_profile.length_params_x[1] = 0.0;
        }
        if (species_param.dens_profile.length_params_x.size()<3) {
            species_param.dens_profile.length_params_x.resize(3);
            species_param.dens_profile.length_params_x[2] = species_param.dens_profile.length_params_x[1];
        }
    }
    
    
    // ------------------------
    // Gaussian density profile
    // ------------------------
    // vacuum_length[0]  : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]    : full length of the density distribution (default value is sim_length-vacuum_length[0])
    // dens_length_x[1]    : FWHM of the gaussian density distribution (default is dens_length_x[0]/3.0)
    // dens_length_x[2]    : center of the gaussian density distribution (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // dens_int_params[0]: order of the gaussian density distribution (default is 2)
    else if (species_param.dens_profile.profile=="gaussian") {
        
        if (species_param.dens_profile.length_params_x.size()<2) {
            species_param.dens_profile.length_params_x.resize(2);
            species_param.dens_profile.length_params_x[1] = species_param.dens_profile.length_params_x[0]/3.0;
        }
        if (species_param.dens_profile.length_params_x.size()<3) {
            species_param.dens_profile.length_params_x.resize(3);
            species_param.dens_profile.length_params_x[2] = species_param.dens_profile.vacuum_length[0]+0.5*species_param.dens_profile.length_params_x[0];
        }
        if (species_param.dens_profile.int_params.size()<1) {
            species_param.dens_profile.int_params.resize(1);
            species_param.dens_profile.int_params[0] = 2;
        }
    }
    
    
    // ---------------------------------------------------------------------------------
    // Charles magnetic field profile for Liang simulations
    // Top-hat profile :
    // int_params[0]      = order of the Gaussian
    // double_params[0]   = background Bfield amplitude
    // double_params[1]   = maximum Bfield amplitude
    // length_params_x[0] = position of the maximum of the B-field
    // length_params_x[1] = FWHM of the magnetic field (Gaussian) distribution
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile == "magexpansion") {
    
        if (species_param.dens_profile.int_params.size()<1)
            ERROR("one int_params must be defined for Charles profile" );
        if (species_param.dens_profile.double_params.size()<2)
            ERROR("two double_params must be defined for Charles profile" );
        if (species_param.dens_profile.length_params_x.size()<2)
            ERROR("three length_params_x must be defined for Charles profile" );
    }
    
    
    // Polygonal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x : specifies the position of different points in the plasma region (0 being the border of the box)
    // dens_dbl_params : specifies the relative (i.e. max=1) value of the density at these points
    else if (species_param.dens_profile.profile=="polygonal") {
        if ( species_param.dens_profile.length_params_x.size() != species_param.dens_profile.double_params.size() )
            ERROR("Incorrect definition of the polygonal profile, size of dens_length_x & dens_dbl_params do not match");
    }
    
    
    // cosine density profile
    // ----------------------
    // vacuum_length[0]   : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the plasma
    // dens_dbl_params[0] : amplitude of the cosine perturbation (has to be <= 1, default = 0.01)
    // dens_dbl_params[1] : number of modes over the plasma length (dens_length_x[0]) (default = 1.0)
    else if (species_param.dens_profile.profile=="cosine") {
        if (species_param.dens_profile.double_params.size()<1) {
            species_param.dens_profile.double_params.resize(1);
            species_param.dens_profile.double_params[0] = 0.01;
        }
        if (species_param.dens_profile.double_params.size()<2) {
            species_param.dens_profile.double_params.resize(2);
            species_param.dens_profile.double_params[1] = 1.0;
        }
        if (species_param.dens_profile.double_params[0]>1.0)
            ERROR("Incorrect definition of the cosine density profile, dens_dbl_params[0] should be <= 1");
    }
    
    
    // Other density profile
    // ---------------------
    else {
        ERROR("Density profile " << species_param.dens_profile.profile << " not defined in 1D");
    }
    
    HEREIAM("");
    
}


double DensityProfile1D::operator() (std::vector<double> x_cell) {
    
    // ------------------------
    // Constant density profile
    // ------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density (default is sim_length-vacuum_length[0])
    if (species_param.dens_profile.profile=="constant") {
        if (   (x_cell[0]>species_param.dens_profile.vacuum_length[0])
            && (x_cell[0]<species_param.dens_profile.vacuum_length[0]+species_param.dens_profile.length_params_x[0]) ) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    
    
    // ---------------------------
    // Trapezoidal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_x[1]   : length of the rising slope (default value is zero)
    // dens_length_x[2]   : length of the descreasing slope (default value is the rising slope value species_length[0])
    else if (species_param.dens_profile.profile=="trapezoidal") {
        
        double vacuum      = species_param.dens_profile.vacuum_length[0];
        double plateau     = species_param.dens_profile.length_params_x[0];
        double left_slope  = species_param.dens_profile.length_params_x[1];
        double right_slope = species_param.dens_profile.length_params_x[2];
 
        // vacuum region
        if ( x_cell[0] < vacuum ) {
            return 0.0;
        }
        // linearly increasing density
        else if ( x_cell[0] < vacuum + left_slope ) {
            return (x_cell[0]-vacuum)/left_slope;
        }
        // density plateau
        else if ( x_cell[0] < vacuum + left_slope + plateau ) {
            return 1.0;
        }
        // linearly decreasing density
        else if ( x_cell[0] < vacuum + left_slope + plateau + right_slope ) {
            return 1.0 - ( x_cell[0] - (vacuum + left_slope + plateau) ) / right_slope;
        }
        // beyond the plasma
        else {
            return 0.0;
        }
    }
    
    
    // ------------------------
    // Gaussian density profile
    // ------------------------
    // vacuum_length[0]  : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]    : full length of the density distribution (default value is sim_length-vacuum_length[0])
    // dens_length_x[1]    : FWHM of the gaussian density distribution (default is dens_length_x[0]/3.0)
    // dens_length_x[2]    : center of the gaussian density distribution (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // dens_int_params[0]: order of the gaussian density distribution (default is 2)
    else if (species_param.dens_profile.profile=="gaussian") {
        
        double vacuum      = species_param.dens_profile.vacuum_length[0];
        double full_length = species_param.dens_profile.length_params_x[0];
        double fwhm        = species_param.dens_profile.length_params_x[1];
        double center      = species_param.dens_profile.length_params_x[2];
        short int N        = species_param.dens_profile.int_params[0];
        double sigmaN      = pow(0.5*fwhm,N)/log(2.0);
        
        // vacuum region
        if ( x_cell[0] < vacuum ) {
            return 0.0;
        }
        // gaussian profile
        else if (x_cell[0] < vacuum+full_length ) {
            return exp( -pow(x_cell[0]-center,N) / sigmaN );
        }
        // beyond the plasma
        else {
            return 0.0;
        }
    }
    
    // ------------------------
    // Charles density profile
    // ------------------------
    // vacuum_length[0]  : not used here
    // int_params[0]      = order of the Gaussian
    // double_params[0]   = background density
    // double_params[1]   = maximum density
    // length_params_x[0] = position of the maximum density
    // length_params_x[1] = FWHM of the density (Gaussian) distribution
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile=="magexpansion") {
        int    N     = species_param.dens_profile.int_params[0];
        double n0    = species_param.dens_profile.double_params[0];
        double nmax  = species_param.dens_profile.double_params[1];
        double x0    = species_param.dens_profile.length_params_x[0];
        double L     = species_param.dens_profile.length_params_x[1];
        double sigma = pow(L/2.0,N)/log(2.0);
        double x     = x_cell[0]-x0;
        
        return (n0 + nmax*exp(-pow(x,N)/sigma))/3.0;
    }
    
    // Polygonal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x : specifies the position of different points in the plasma region (0 being the border of the box)
    // dens_dbl_params : specifies the relative (i.e. max=1) value of the density at these points
    else if (species_param.dens_profile.profile=="polygonal") {
        
        unsigned int N = species_param.dens_profile.length_params_x.size();
        
        // vacuum region
        if ( x_cell[0] < species_param.dens_profile.vacuum_length[0] ) {
            return 0.0;
        }
        // plasma region (defined over N segments)
        else if(x_cell[0] < species_param.dens_profile.length_params_x[N] ){
            for( unsigned int i=1; i<N; ++i){
                if ( (x_cell[0]>=species_param.dens_profile.length_params_x[i-1]) && (x_cell[0]<species_param.dens_profile.length_params_x[i]) ){
                    double m = (species_param.dens_profile.double_params[i]-species_param.dens_profile.double_params[i-1])
                    /          (species_param.dens_profile.length_params_x[i]-species_param.dens_profile.length_params_x[i-1]);
                    
                    return species_param.dens_profile.double_params[i-1] + m * ( x_cell[0]-species_param.dens_profile.length_params_x[i-1] );
                }
            }
        }
        // beyond the plasma
        else {
            return 0.0;
        }
    }
    
    
    // cosine density profile
    // ----------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the plasma
    // dens_dbl_params[0]   : amplitude of the cosine perturbation (has to be <= 1)
    // dens_dbl_params[1]   : number of modes over the plasma length (dens_length_x[0])
    else if (species_param.dens_profile.profile=="cosine") {
        
        // vacuum region
        if ( x_cell[0] < species_param.dens_profile.vacuum_length[0] ) {
            return 0.0;
        }
        // plasma region
        else if (x_cell[0] < species_param.dens_profile.vacuum_length[0]+species_param.dens_profile.length_params_x[0]) {
            double x = x_cell[0]-species_param.dens_profile.vacuum_length[0];
            double k = 2.0*M_PI * species_param.dens_profile.double_params[1] / species_param.dens_profile.length_params_x[0];
            return 1.0 + species_param.dens_profile.double_params[0] * cos(k*x);
        }
        // beyond the plasma
        else {
            return 0.0;
        }
    }
    
    
    // Other density profile
    // ---------------------
    else {
        ERROR("Density profile " << species_param.dens_profile.profile << " not defined in 1D");
    }
    
    return 0;
};
