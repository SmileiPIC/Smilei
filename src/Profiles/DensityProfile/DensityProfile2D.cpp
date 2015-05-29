#include "DensityProfile2D.h"
#include "Tools.h"

using namespace std;

DensityProfile2D::DensityProfile2D(SpeciesStructure &params) : DensityProfile(params) {
    
    // Constant density profile
    // ------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
    // species_length_x[0]: length of the plasma in the x-direction
    // species_length_y[0]: length of the plasma in the y-direction
    if (species_param.dens_profile.profile=="constant") {
        // nothing to be done here, all default parameters are computed directly in PicParams.cpp
    }
    
    // Trapezoidal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_x[1]   : length of the left slope (default value is zero)
    // dens_length_x[2]   : length of the right slope (default value is the rising slope value species_length[1])
    // same definition of dens_length_y (but in y-direction)
    else if (species_param.dens_profile.profile=="trapezoidal") {
        
        // x-direction
        if (species_param.dens_profile.length_params_x.size()<2) {
            species_param.dens_profile.length_params_x.resize(2);
            species_param.dens_profile.length_params_x[1] = 0.0;
        }
        if (species_param.dens_profile.length_params_x.size()<3) {
            species_param.dens_profile.length_params_x.resize(3);
            species_param.dens_profile.length_params_x[2] = species_param.dens_profile.length_params_y[1];
        }
        //y-direction
        if (species_param.dens_profile.length_params_y.size()<2) {
            species_param.dens_profile.length_params_y.resize(2);
            species_param.dens_profile.length_params_y[1] = 0.0;
        }
        if (species_param.dens_profile.length_params_y.size()<3) {
            species_param.dens_profile.length_params_y.resize(3);
            species_param.dens_profile.length_params_y[2] = species_param.dens_profile.length_params_y[1];
        }
        
    }
    
    // Gaussian density profile
    // ------------------------
    // vacuum_length[0]  : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]  : full length of the density distribution (default value is sim_length-vacuum_length[0])
    // dens_length_x[1]  : FWHM of the gaussian density distribution (default is dens_length_x[0]/3.0)
    // dens_length_x[2]  : center of the gaussian density distribution (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // dens_int_params[0]: order of the gaussian density distribution (default is 2)
    // same definitions hold for the y-direction with dens_int_params[1] the order of the Gaussian
    // note that if dens_int_params[1]=0 (default value) then the profile is constant in the y-direction
    else if (species_param.dens_profile.profile=="gaussian") {
        
        // x-direction
        if (species_param.dens_profile.int_params.size()<1) {
            species_param.dens_profile.int_params.resize(1);
            species_param.dens_profile.int_params[0] = 2;
        }
        if (species_param.dens_profile.length_params_x.size()<2) {
            species_param.dens_profile.length_params_x.resize(2);
            species_param.dens_profile.length_params_x[1] = species_param.dens_profile.length_params_x[0]/3.0;
        }
        if (species_param.dens_profile.length_params_x.size()<3) {
            species_param.dens_profile.length_params_x.resize(3);
            species_param.dens_profile.length_params_x[2] = species_param.dens_profile.vacuum_length[0]+0.5*species_param.dens_profile.length_params_x[0];
        }
        
        // y-direction
        if (species_param.dens_profile.int_params.size()<2) {
            species_param.dens_profile.int_params.resize(2);
            species_param.dens_profile.int_params[0] = 0;
        }
        if (species_param.dens_profile.length_params_y.size()<2) {
            species_param.dens_profile.length_params_y.resize(2);
            species_param.dens_profile.length_params_y[1] = species_param.dens_profile.length_params_y[0]/3.0;
        }
        if (species_param.dens_profile.length_params_y.size()<3) {
            species_param.dens_profile.length_params_y.resize(3);
            species_param.dens_profile.length_params_y[2] = species_param.dens_profile.vacuum_length[1]+0.5*species_param.dens_profile.length_params_y[0];
        }
        
    }
    
    else if (species_param.dens_profile.profile=="preplasma") {
        // Preplasma density profile (exponential pre/post plasma in x-direction, constant in y)
        // -------------------------------------------------------------------------------------
        // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
        // species_length_x[0]: length of the preplasma
        // species_length_x[1]: slope of the preplasma
        // species_length_x[2]: length of the plateau
        // species_length_x[3]: length of the postplasma
        // species_length_x[4]: slope of the postplasma
        // species_length_y[0]: length of the plasma in the y-direction
        if (species_param.dens_profile.length_params_x.size()<5) {
            ERROR("For the Preplasma density profile 5 length_params_x have to be defined");
        }

    }

    
    // Harris density profile: used for reconnection
    // ---------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[0] : characteristic length of the Harris profile
    // double_params[0]   : density nb parameter
    else if (species_param.dens_profile.profile=="harris") {
        
        if (species_param.dens_profile.double_params.size()<1) {
            ERROR("For the Harris density profile 1 double_params has to be defined");
        }
        if (species_param.dens_profile.length_params_y.size()<3) {
            ERROR("For the Harris density profile 3 length_params_y have to be defined");
        }
        
    }
    
    // ---------------------------------------------------------------------------------
    // Charles magnetic field profile for Liang simulations
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_y[0] = position of the maximum magnetic field
    // length_params_y[1] = Length of the magnetic gradient
    // dens_length_y[2]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_y[3]   : length of the left slope (default value is zero)
    // dens_length_y[4]   : length of the right slope (default value is the rising slope value species_length[1])
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile == "magexpansion") {
    
        //if (species_param.dens_profile.int_params.size()<1)
        //    ERROR("one int_params must be defined for Charles profile" );
        if (species_param.dens_profile.double_params.size()<4)
            ERROR("two double_params must be defined for Charles profile" );
        if (species_param.dens_profile.length_params_y.size()<2)
            ERROR("two length_params_y must be defined for Charles profile" );
	    
	 
        //y-direction
        if (species_param.dens_profile.length_params_y.size()<3) {
            species_param.dens_profile.length_params_y.resize(3);
            species_param.dens_profile.length_params_y[2] = 1.e+10;
        }
        if (species_param.dens_profile.length_params_y.size()<4) {
            species_param.dens_profile.length_params_y.resize(4);
            species_param.dens_profile.length_params_y[3] = 0.0;
        }
        if (species_param.dens_profile.length_params_y.size()<5) {
            species_param.dens_profile.length_params_y.resize(5);
            species_param.dens_profile.length_params_y[4] = species_param.dens_profile.length_params_y[3];
        }
    }
    
    // ---------------------------------------------------------------------------------
    // Blob magnetic field profile for Korneev simulations
    // vacuum_length[0]  : not used here
    // double_params[0] = Relative Density variation 
    // dens_length_x[0] = x position of the maximum magnetic field
    // dens_length_x[1] = Length of the density gradient
    // dens_length_y[0] = y position of the maximum magnetic field
    // dens_length_y[1]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_y[2]   : length of the left slope (default value is zero)
    // dens_length_y[3]   : length of the right slope (default value is the rising slope value species_length[1])
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile == "blob") {
    
        if (species_param.dens_profile.double_params.size()<1)
            ERROR("two double_params must be defined for Charles profile" );
        if (species_param.dens_profile.length_params_y.size()<1)
            ERROR("two length_params_y must be defined for Charles profile" );
        if (species_param.dens_profile.length_params_x.size()<2)
            ERROR("two length_params_x must be defined for Charles profile" );
	   
	 
        //y-direction
        if (species_param.dens_profile.length_params_y.size()<2) {
            species_param.dens_profile.length_params_y.resize(2);
            species_param.dens_profile.length_params_y[1] = 1.e+10;
        }
        if (species_param.dens_profile.length_params_y.size()<3) {
            species_param.dens_profile.length_params_y.resize(3);
            species_param.dens_profile.length_params_y[2] = 0.0;
        }
        if (species_param.dens_profile.length_params_y.size()<4) {
            species_param.dens_profile.length_params_y.resize(4);
            species_param.dens_profile.length_params_y[3] = species_param.dens_profile.length_params_y[2];
        }
    }
    
    // Grating
    // -------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // double_params[0]   : phase of the cosine perturbation at the center of the grating
    // length_params_x[0] : full length of the plasma in the x-direction
    // length_params_x[1] : depth of the grating
    // length_params_y[0] : full length of the plasma in the y-direction
    // length_params_y[1] : period of the grating
    else if (species_param.dens_profile.profile=="grating") {
        
        if (species_param.dens_profile.double_params.size()<1) {
            species_param.dens_profile.double_params.resize(1);
            species_param.dens_profile.double_params[0] = 0.0;
            WARNING("For the Grating density profile phase double_params[0] put to zero");
        }
        if (   (species_param.dens_profile.length_params_x.size()<2)
            || (species_param.dens_profile.length_params_y.size()<2) ) {
            ERROR("For the Grating density profile at least 2 parameters length_params_x/y must be defined");
        }
        
    }    
    else if (species_param.dens_profile.profile=="python") {
        DEBUG("it's a python profile");
    }
    else {
        ERROR("Profile " << species_param.dens_profile.profile << " not defined");
    }//if profile

    
}

double DensityProfile2D::operator() (vector<double> x_cell) {
    
    double fx, fy;
    
    // Constant density profile
    // ------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
    // species_length_x[0]: length of the plasma in the x-direction
    // species_length_y[0]: length of the plasma in the y-direction
    if (species_param.dens_profile.profile=="constant") {
        
        // x-direction
        if (   (x_cell[0]>species_param.dens_profile.vacuum_length[0])
            && (x_cell[0]<species_param.dens_profile.vacuum_length[0]+species_param.dens_profile.length_params_x[0]) ) {
            fx = 1.0;
        } else {
            fx = 0.0;
        }
        
        // y-direction
        if (   (x_cell[1]>species_param.dens_profile.vacuum_length[1])
            && (x_cell[1]<species_param.dens_profile.vacuum_length[1]+species_param.dens_profile.length_params_y[0]) ) {
            fy = 1.0;
        } else {
            fy = 0.0;
        }
        
        // x-y direction
        return fx*fy;
        
    }// constant
    
    
    // Trapezoidal density profile
    // ---------------------------
    // vacuum_length[0] : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_x[1]   : length of the left slope (default value is zero)
    // dens_length_x[2]   : length of the right slope (default value is the rising slope value species_length[1])
    // same definition of dens_length_y (but in y-direction)
    else if (species_param.dens_profile.profile=="trapezoidal") {
        
        // x-direction
        double vacuum      = species_param.dens_profile.vacuum_length[0];
        double plateau     = species_param.dens_profile.length_params_x[0];
        double left_slope  = species_param.dens_profile.length_params_x[1];
        double right_slope = species_param.dens_profile.length_params_x[2];
        
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
        
        // x-direction
        vacuum      = species_param.dens_profile.vacuum_length[1];
        plateau     = species_param.dens_profile.length_params_y[0];
        left_slope  = species_param.dens_profile.length_params_y[1];
        right_slope = species_param.dens_profile.length_params_y[2];
        
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
        
        // x-y directions
        return fx*fy;
        
    }// trapezoidal
    
    // ------------------------
    // Charles density profile
    // ------------------------
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_x[0] = position of the maximum magnetic field
    // length_params_x[1] = Length of the magnetic gradient
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile=="magexpansion") {
        //int    N     = species_param.dens_profile.int_params[0];
        //double n0    = species_param.dens_profile.double_params[0];
        double P0    = species_param.dens_profile.double_params[1];
        double B0    = species_param.dens_profile.double_params[2];
        double Bmax  = species_param.dens_profile.double_params[3];
        double x0    = species_param.dens_profile.length_params_y[0];
        double L     = species_param.dens_profile.length_params_y[1];
	double tiny  = 1e-10*L;
        //double sigma = pow(L/2.0,N)/log(2.0);
        double x     = x_cell[1]-x0;
	
	double fy;
        
        // y-direction
        double vacuum      = species_param.dens_profile.vacuum_length[1];
        double plateau     = species_param.dens_profile.length_params_y[2];
        double left_slope  = species_param.dens_profile.length_params_y[3];
        double right_slope = species_param.dens_profile.length_params_y[4];
        
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
	
	
	
	if (Bmax == 0.) {
		double Bm = sqrt(pow(B0,2) + 2*P0)-B0;
		double B  = B0 + Bm/pow(cosh(x/L),2);
		double A  = B0*x + Bm*L*tanh(x/L);
		double DP = P0 + pow(B0,2)/2 - pow(B,2)/2;
		if (abs(x)<tiny) {
			return fy*(exp(-2));
			//return 1.;
		}
		else {
        		return fy*(exp( -2*A*Bm/L*tanh(x/L) /(DP*pow(cosh(x/L),2)) ));
        		//return (exp( - abs(tanh(x/L) )));
			//return 1.;
		}
	}}
	
    // ---------------------------------------------------------------------------------
    // Blob magnetic field profile for Korneev simulations
    // vacuum_length[0]  : used here
    // double_params[0]   = Density variation 
    // dens_length_x[0] = x position of the maximum magnetic field
    // dens_length_x[1] = Length of the density gradient
    // dens_length_y[0] = y position of the maximum magnetic field
    // dens_length_y[1]   : length of the density plateau (default value if sim_length-vacuum_length[0] as for constant)
    // dens_length_y[2]   : length of the left slope (default value is zero)
    // dens_length_y[3]   : length of the right slope (default value is the rising slope value species_length[1])
    // ---------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile=="blob") {
        double dn    = species_param.dens_profile.double_params[0];
        double x0    = species_param.dens_profile.length_params_x[0];
        double y0    = species_param.dens_profile.length_params_y[0];
        double L     = species_param.dens_profile.length_params_x[1];
        double r     = sqrt(pow(x_cell[0]-x0,2) + pow(x_cell[1]-y0,2));
	
	double fy;
        
        // y-direction
        double vacuum      = species_param.dens_profile.vacuum_length[1];
        double plateau     = species_param.dens_profile.length_params_y[1];
        double left_slope  = species_param.dens_profile.length_params_y[2];
        double right_slope = species_param.dens_profile.length_params_y[3];
        
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
	
	
        return fy*( 1 - dn/(pow(cosh(r/L),2)) );
	}
    // Gaussian profile
    // ----------------
    // vacuum_length[0]  : length of the vacuum region before the plasma (default is 0)
    // dens_length_x[0]  : full length of the density distribution (default value is sim_length-vacuum_length[0])
    // dens_length_x[1]  : FWHM of the gaussian density distribution (default is dens_length_x[0]/3.0)
    // dens_length_x[2]  : center of the gaussian density distribution (where it is maximum)
    //                     (default is vaccum_length + 1/2 of full length)
    // dens_int_params[0]: order of the gaussian density distribution (default is 2)
    // same definitions hold for the y-direction with dens_int_params[1] the order of the Gaussian
    // note that if dens_int_params[1]=0 (default value) then the profile is constant in the y-direction
    else if (species_param.dens_profile.profile=="gaussian") {
        
        // x-direction
        short int N        = species_param.dens_profile.int_params[0];
        double vacuum      = species_param.dens_profile.vacuum_length[0];
        double full_length = species_param.dens_profile.length_params_x[0];
        double fwhm        = species_param.dens_profile.length_params_x[1];
        double center      = species_param.dens_profile.length_params_x[2];
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
        N  = species_param.dens_profile.int_params[1];
        fy = 1.0;
        
        if (N>0) { // if N=0, then returns constant profile in the y-direction (fy=1.0)
            vacuum      = species_param.dens_profile.vacuum_length[1];
            full_length = species_param.dens_profile.length_params_y[0];
            fwhm        = species_param.dens_profile.length_params_y[1];
            center      = species_param.dens_profile.length_params_y[2];
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
        
        // x-y directions
        return fx*fy;
        
    }// gaussian
    
    // Preplasma density profile (exponential pre/post plasma in x-direction, constant in y)
    // -------------------------------------------------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
    // species_length_x[0]: length of the preplasma
    // species_length_x[1]: slope of the preplasma
    // species_length_x[2]: length of the plateau
    // species_length_x[3]: length of the postplasma
    // species_length_x[4]: slope of the postplasma
    // species_length_y[0]: length of the plasma in the y-direction
    if (species_param.dens_profile.profile=="preplasma") {
        
        double l0 = species_param.dens_profile.vacuum_length[0];
        double lf = species_param.dens_profile.length_params_x[0];
        double l1 = species_param.dens_profile.length_params_x[1];
        double lp = species_param.dens_profile.length_params_x[2];
        double lb = species_param.dens_profile.length_params_x[3];
        double l2 = species_param.dens_profile.length_params_x[4];
        
        // x-direction
        if (x_cell[0]<l0) {
	    fx = 0.0;
        } else if (x_cell[0]<l0+lf) {
            fx = exp( (x_cell[0]-(l0+lf))/l1 );
        } else if (x_cell[0]<l0+lf+lp) {
            fx = 1.0;
        } else if (x_cell[0]<l0+lf+lp+lb) {
            fx = exp( -(x_cell[0]-(l0+lf+lp))/l2 );
        } else {
            fx = 0.0;
        }
        
        // y-direction
        if (   (x_cell[1]>species_param.dens_profile.vacuum_length[1])
            && (x_cell[1]<species_param.dens_profile.vacuum_length[1]+species_param.dens_profile.length_params_y[0]) ) {
            fy = 1.0;
        } else {
            fy = 0.0;
        }
        
        // x-y direction
        return fx*fy;
        
    }// constant

    // Harris density profile: used for reconnection
    // ---------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[1] : position of the first density maximum
    // length_params_y[2] : position of the second density maximum
    // double_params[0]   : density nb parameter
    else if (species_param.dens_profile.profile=="harris") {
    
        double nb = species_param.dens_profile.double_params[0];
        double L  = species_param.dens_profile.length_params_y[0];
        double y0 = species_param.dens_profile.length_params_y[1];
        double y1 = species_param.dens_profile.length_params_y[2];
        
        return 1.0 + 1.0/( nb * pow(cosh((x_cell[1]-y0)/L),2) ) + 1.0/( nb * pow(cosh((x_cell[1]-y1)/L),2) );
        
    }
    
    // Grating
    // -------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // double_params[0]   : phase of the cosine perturbation at the center of the grating
    // length_params_x[0] : full length of the plasma in the x-direction
    // length_params_x[1] : depth of the grating
    // length_params_y[0] : full length of the plasma in the y-direction
    // length_params_y[1] : period of the grating
    else if (species_param.dens_profile.profile=="grating") {
        
        double x0    = species_param.dens_profile.vacuum_length[0];
        double y0    = species_param.dens_profile.vacuum_length[1];
        double Lx    = species_param.dens_profile.length_params_x[0];
        double delta = species_param.dens_profile.length_params_x[1];
        double Ly    = species_param.dens_profile.length_params_y[0];
        double k     = 2.0*M_PI / species_param.dens_profile.length_params_y[1];
        double phi   = species_param.dens_profile.double_params[0];
        // position of the border
        double xb = x0 + delta * cos( k*(x_cell[1]-(y0+Ly/2.0)) + phi );
        
        if (   (xb<x_cell[0]) && (x_cell[0]<x0+Lx)
            && (y0<x_cell[1]) && (x_cell[1]<y0+Ly) ) {
            return 1.0;
        } else {
            return 0.0;
        }
        
    }
    
    // Plasma density profile corresponding to Fukuda et al., Phys. Rev. Lett. 103, 165002 (2012)
    // used in simulations for Anna Levy
    // ------------------------------------------------------------------------------------------
    else if (species_param.dens_profile.profile=="fukuda") {
        
        // x-direction
        if (x_cell[0]<2.0*M_PI*2.0) {
            fx = 0.0;
        }
        else if (x_cell[0]<2.0*M_PI*13.0) {
            fx = 0.2;
        }
        else if (x_cell[0]<2.0*M_PI*20.0) {
            fx = 0.2 + 0.8*(x_cell[0]-2.0*M_PI*13.0)/(2.0*M_PI*7.0);
        }
        else if (x_cell[0]<2.0*M_PI*65.0) {
            fx = 1.0;
        }
        else if (x_cell[0]<2.0*M_PI*82.0) {
            fx = 1.0 - 0.8*(x_cell[0]-2.0*M_PI*65.0)/(2.0*M_PI*17.0);
        }
        else if (x_cell[0]<2.0*M_PI*112.0) {
            fx = 0.2;
        }
        else {
            fx = 0.0;
        }
        
        // y-direction: constant density
        fy = 1.0;
        
        // x-y direction
        return fx*fy;
        
    }// fukuda
    else if (species_param.dens_profile.profile=="python") {
        return PyTools::runPyFunction(species_param.dens_profile.py_profile, x_cell[0], x_cell[1]);
    }
    // Other profiles: not defined
    // ---------------------------
    else {
        ERROR("Density profile '" << species_param.dens_profile.profile <<"' not yet defined in 2D");
        return 0.0;
    }

}

