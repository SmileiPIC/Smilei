#include "Profile.h"


using namespace std;

// Default constructor.
// Applies to profiles for species (density, velocity and temperature profiles)
Profile::Profile(ProfileStructure & pp, string geometry, double convfac) :
conv_fac(convfac)
{
    
    factor = 1.;
    profile_param = pp;
    
    // Launch the initialization
    init(profile_param, geometry);
    
}


// Special constructor.
// Applies to external field profiles
Profile::Profile(ExtFieldStructure & pp, string geometry, double convfac) :
conv_fac(convfac)
{
    profile_param = static_cast<ProfileStructure> (pp);
    
    // define the factor
    factor = pp.factor;
    // define the vacuum_length as zeros
    profile_param.vacuum_length.resize(dim);
    for( int i=0; i<dim; i++) profile_param.vacuum_length[i]=0.;
    
    // Launch the initialization
    init(profile_param, geometry);
    
}


void Profile::init(ProfileStructure & pp, string geometry)
{

    if      (geometry == "1d3v") dim = 1;
    else if (geometry == "2d3v") dim = 2;
    else {
        ERROR( "Unsupported geometry : " << geometry);
    }
    
    
    // For all profiles:
    // vacuum_length : length of the vacuum region before the start of the profile (default is 0)
    //                 Vector with 1, 2 or 3, components, corresponding to X, Y, Z.
    
    // ------------------------
    // Constant profile
    // ------------------------
    // ****_length_x[0] : length of the profile (default is sim_length-vacuum_length[0])
    if (profile_param.profile=="constant") {
        // nothing done here, by default: vacuum_length[0] = 0, ****_length_x[0] = 0
    }
    
    // ---------------------------
    // Trapezoidal profile
    // ---------------------------
    // ****_length_x[0] : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // ****_length_x[1] : length of the left slope (default value is zero)
    // ****_length_x[2] : length of the right slope (default value is the rising slope value species_length[1])
    // same definition of ****_length_y (but in y-direction)
    else if (profile_param.profile=="trapezoidal") {
        // X
        if (profile_param.length_params_x.size()<2) {
            profile_param.length_params_x.resize(2);
            profile_param.length_params_x[1] = 0.0;
        }
        if (profile_param.length_params_x.size()<3) {
            profile_param.length_params_x.resize(3);
            profile_param.length_params_x[2] = profile_param.length_params_x[1];
        }
        
        // Y
        if( dim >= 2 ) {
            if (profile_param.length_params_y.size()<2) {
                profile_param.length_params_y.resize(2);
                profile_param.length_params_y[1] = 0.0;
            }
            if (profile_param.length_params_y.size()<3) {
                profile_param.length_params_y.resize(3);
                profile_param.length_params_y[2] = profile_param.length_params_y[1];
            }
        }
    }
    
    
    // ------------------------
    // Gaussian profile
    // ------------------------
    // ****_length_x[0]  : full length of the distribution (default value is sim_length-vacuum_length[0])
    // ****_length_x[1]  : FWHM of the gaussian (default is ****_length_x[0]/3.0)
    // ****_length_x[2]  : center of the gaussian (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // ****_int_params[0]: order of the gaussian (default is 2)
    // same definitions hold for the y-direction with ***_int_params[1] the order of the Gaussian
    // Note that if ****_int_params[1]=0 (default value) then the profile is constant in the y-direction
    // same definition of ****_length_y (but in y-direction)
    else if (profile_param.profile=="gaussian") {
        // X
        if (profile_param.length_params_x.size()<2) {
            profile_param.length_params_x.resize(2);
            profile_param.length_params_x[1] = profile_param.length_params_x[0]/3.0;
        }
        if (profile_param.length_params_x.size()<3) {
            profile_param.length_params_x.resize(3);
            profile_param.length_params_x[2] = profile_param.vacuum_length[0]+0.5*profile_param.length_params_x[0];
        }
        if (profile_param.int_params.size()<1) {
            profile_param.int_params.resize(1);
            profile_param.int_params[0] = 2;
        }
        
        // Y
        if( dim >= 2 ) {
            if (profile_param.int_params.size()<2) {
                profile_param.int_params.resize(2);
                profile_param.int_params[0] = 0;
            }
            if (profile_param.length_params_y.size()<2) {
                profile_param.length_params_y.resize(2);
                profile_param.length_params_y[1] = profile_param.length_params_y[0]/3.0;
            }
            if (profile_param.length_params_y.size()<3) {
                profile_param.length_params_y.resize(3);
                profile_param.length_params_y[2] = profile_param.vacuum_length[1]+0.5*profile_param.length_params_y[0];
            }
        }
    }
    
    
    // ---------------------------
    // Polygonal profile
    // ---------------------------
    // ****_length_x    : specifies the position of different points of the profile (0 being the border of the box)
    // ****_dbl_params  : specifies the relative (i.e. max=1) value of the profile at these points
    // Note: Polygon applies only to X
    else if (profile_param.profile=="polygonal") {
        if ( profile_param.length_params_x.size() != profile_param.double_params.size() ) {
            ERROR("Incorrect definition of the polygonal profile, size of ****_length_x & ****_dbl_params do not match");
        }
    }
    
    
    // ---------------------------
    // cosine profile
    // ----------------------
    // ****_length_x[0]   : length of the profile
    // ****_dbl_params[0] : amplitude of the cosine perturbation (has to be <= 1, default = 0.01)
    // ****_dbl_params[1] : number of periods over the profile length (****_length_x[0]) (default = 1.0)
    // Note: Cosine applies only to X
    else if (profile_param.profile=="cosine") {
        if (profile_param.double_params.size()<1) {
            profile_param.double_params.resize(1);
            profile_param.double_params[0] = 0.01;
        }
        if (profile_param.double_params.size()<2) {
            profile_param.double_params.resize(2);
            profile_param.double_params[1] = 1.0;
        }
        if (profile_param.double_params[0]>1.0) {
            ERROR("Incorrect definition of the cosine density profile, ****_dbl_params[0] should be <= 1");
        }
    }
    
    
    else if (profile_param.profile=="python") {
        DEBUG("it's a python profile");
    }
    
    
    else {
        ERROR("Profile " << profile_param.profile << " not defined");
    }

}


Profile::~Profile()
{
}


double Profile::valueAt (vector<double> x_cell) {
    
    double fx=1., fy=1.;
    
    // ------------------------
    // Constant profile
    // ------------------------
    // ****_length_x[0] : length of the profile (default is sim_length-vacuum_length[0])
    if (profile_param.profile=="constant") {
        // x-direction
        if (   (x_cell[0]<profile_param.vacuum_length[0])
            || (x_cell[0]>profile_param.vacuum_length[0]+profile_param.length_params_x[0]) ) {
            fx = 0.;
        }
        // y-direction
        if( dim >= 2 ) {
            if (   (x_cell[1]<profile_param.vacuum_length[1])
                || (x_cell[1]>profile_param.vacuum_length[1]+profile_param.length_params_y[0]) ) {
                fy = 0.;
            }
        }
        
        // x-y direction
        return factor*fx*fy;
    }
    
    
    // ---------------------------
    // Trapezoidal profile
    // ---------------------------
    // ****_length_x[0] : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // ****_length_x[1] : length of the left slope (default value is zero)
    // ****_length_x[2] : length of the right slope (default value is the rising slope value species_length[1])
    // same definition of ****_length_y (but in y-direction)
    else if (profile_param.profile=="trapezoidal") {
        
        // x-direction
        double vacuum      = profile_param.vacuum_length[0];
        double plateau     = profile_param.length_params_x[0];
        double left_slope  = profile_param.length_params_x[1];
        double right_slope = profile_param.length_params_x[2];
        
        // vacuum region
        if ( x_cell[0] < vacuum ) {
            fx = 0.0;
        }
        // linearly increasing density
        else if ( x_cell[0] < vacuum+left_slope ) {
            fx = (x_cell[0]-vacuum) / left_slope;
        }
        // density plateau
        else if ( x_cell[0] < vacuum+left_slope+plateau ) {
            fx = 1.0;
        }
        // linearly decreasing density
        else if ( x_cell[0] < vacuum+left_slope+plateau+right_slope ) {
            fx = 1.0 - ( x_cell[0] - (vacuum+left_slope+right_slope) ) / right_slope;
        }
        // beyond the plasma
        else {
            fx = 0.0;
        }
        
        // y-direction
        if( dim >= 2 ) {
            vacuum      = profile_param.vacuum_length[1];
            plateau     = profile_param.length_params_y[0];
            left_slope  = profile_param.length_params_y[1];
            right_slope = profile_param.length_params_y[2];
            
            // vacuum region
            if ( x_cell[1] < vacuum ) {
                fy = 0.0;
            }
            // linearly increasing density
            else if ( x_cell[1] < vacuum+left_slope ) {
                fy = (x_cell[1]-vacuum) / left_slope;
            }
            // density plateau
            else if ( x_cell[1] < vacuum+left_slope+plateau ) {
                fy = 1.0;
            }
            // linearly decreasing density
            else if ( x_cell[1] < vacuum+left_slope+plateau+right_slope ) {
                fy = 1.0 - ( x_cell[1] - (vacuum+left_slope+right_slope) ) / right_slope;
            }
            // beyond the plasma
            else {
                fy = 0.0;
            }
        }
        
        // x-y directions
        return factor*fx*fy;
    }
    
    
    // ------------------------
    // Gaussian profile
    // ------------------------
    // ****_length_x[0]  : full length of the distribution (default value is sim_length-vacuum_length[0])
    // ****_length_x[1]  : FWHM of the gaussian (default is ****_length_x[0]/3.0)
    // ****_length_x[2]  : center of the gaussian (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // ****_int_params[0]: order of the gaussian (default is 2)
    // same definitions hold for the y-direction with ***_int_params[1] the order of the Gaussian
    // Note that if ****_int_params[1]=0 (default value) then the profile is constant in the y-direction
    // same definition of ****_length_y (but in y-direction)
    else if (profile_param.profile=="gaussian") {
        
        // x-direction
        short int N        = profile_param.int_params[0];
        double vacuum      = profile_param.vacuum_length[0];
        double full_length = profile_param.length_params_x[0];
        double fwhm        = profile_param.length_params_x[1];
        double center      = profile_param.length_params_x[2];
        double sigmaN      = pow(0.5*fwhm,N)/log(2.0);
        
        // vacuum region
        if ( x_cell[0] < vacuum ) {
            fx = 0.0;
        }
        // gaussian profile
        else if (x_cell[0] < vacuum+full_length ) {
            fx = exp( -pow(x_cell[0]-center,N) / sigmaN );
        }
        // beyond the plasma
        else {
            fx = 0.0;
        }
        
        // y-direction
        if( dim >= 2 ) {
            N  = profile_param.int_params[1];
            if (N>0) { // if N=0, then returns constant profile in the y-direction (fy=1.0)
                vacuum      = profile_param.vacuum_length[1];
                full_length = profile_param.length_params_y[0];
                fwhm        = profile_param.length_params_y[1];
                center      = profile_param.length_params_y[2];
                sigmaN      = pow(0.5*fwhm,N)/log(2.0);
                
                // vacuum region
                if ( x_cell[1] < vacuum ) {
                    fy = 0.0;
                }
                // gaussian profile
                else if (x_cell[1] < vacuum+full_length ) {
                    fy = exp( -pow(x_cell[1]-center,N) / sigmaN );
                }
                // beyond the plasma
                else {
                    fy = 0.0;
                }
            }
        }
        
        // x-y directions
        return factor*fx*fy;
    }
    
    
    // ---------------------------
    // Polygonal profile
    // ---------------------------
    // ****_length_x    : specifies the position of different points of the profile (0 being the border of the box)
    // ****_dbl_params  : specifies the relative (i.e. max=1) value of the profile at these points
    // Note: Polygon applies only to X
    else if (profile_param.profile=="polygonal") {
        
        unsigned int N = profile_param.length_params_x.size();
        
        // vacuum region
        if ( x_cell[0] < profile_param.vacuum_length[0] ) {
            return 0.0;
        }
        // plasma region (defined over N segments)
        else if(x_cell[0] < profile_param.length_params_x[N] ){
            for( unsigned int i=1; i<N; ++i){
                if ( (x_cell[0]>=profile_param.length_params_x[i-1]) && (x_cell[0]<profile_param.length_params_x[i]) ){
                    double m = (profile_param.double_params  [i]-profile_param.double_params  [i-1])
                              /(profile_param.length_params_x[i]-profile_param.length_params_x[i-1]);
                    
                    return factor*(profile_param.double_params[i-1] + m * ( x_cell[0]-profile_param.length_params_x[i-1] ));
                }
            }
        }
        // beyond the plasma
        else {
            return 0.0;
        }
    }
    
    
    // ----------------------
    // cosine profile
    // ----------------------
    // ****_length_x[0]   : length of the profile
    // ****_dbl_params[0] : amplitude of the cosine perturbation (has to be <= 1, default = 0.01)
    // ****_dbl_params[1] : number of periods over the profile length (****_length_x[0]) (default = 1.0)
    // Note: Cosine applies only to X
    else if (profile_param.profile=="cosine") {
        
        // vacuum region
        if ( x_cell[0] < profile_param.vacuum_length[0] ) {
            return 0.0;
        }
        // profile region
        else if (x_cell[0] < profile_param.vacuum_length[0]+profile_param.length_params_x[0]) {
            double x = x_cell[0]-profile_param.vacuum_length[0];
            double k = 2.0*M_PI * profile_param.double_params[1] / profile_param.length_params_x[0];
            return factor*(1.0 + profile_param.double_params[0] * cos(k*x));
        }
        // beyond the profile
        else {
            return 0.0;
        }
    }
    
    
    else if (profile_param.profile=="python") {
        if        ( dim == 1 ) {
            return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]/conv_fac);
        } else if ( dim == 2 ) {
            return PyTools::runPyFunction(profile_param.py_profile, x_cell[0]/conv_fac, x_cell[1]/conv_fac);
        }
    }
    
    // Other density profile
    // ---------------------
    else {
        ERROR("Density profile " << profile_param.profile << " not defined in 1D");
    }
    
    return 0;
};
