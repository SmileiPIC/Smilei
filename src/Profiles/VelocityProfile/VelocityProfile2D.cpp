#include "VelocityProfile2D.h"
#include "Tools.h"

using namespace std;

VelocityProfile2D::VelocityProfile2D(ProfileSpecies &my_prof_params) : VelocityProfile(my_prof_params) {
    
    // Constant velocity profile
    // -------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
    // species_length_x[0]: length of the plasma in the x-direction
    // species_length_y[0]: length of the plasma in the y-direction
    if (prof_params.profile=="constant") {
        // nothing to be done here, all default parameters are computed directly in PicParams.cpp
    }
    
    
    // Harris velocity profile: used for reconnection
    // ----------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[0] : characteristic length of the Harris profile
    // double_params[0]   : density nb parameter
    else if (prof_params.profile=="harris") {
        
        if (prof_params.double_params.size()<1) {
            ERROR("For the Harris velocity profile 1 double_params has to be defined");
        }
        if (prof_params.length_params_y.size()<2) {
            ERROR("For the Harris velocity profile 2 length_params_y have to be defined");
        }
    }
    
    else {
        ERROR("Profile " << prof_params.profile << " is not defined");
    }//if species_geometry
    
}

double VelocityProfile2D::operator() (vector<double> x_cell) {
    double fx, fy;
    
    // Constant velocity profile
    // -------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions, respectively
    // species_length_x[0]: length of the plasma in the x-direction
    // species_length_y[0]: length of the plasma in the y-direction
    if (prof_params.profile=="constant") {
        
        // x-direction
        if (   (x_cell[0]>prof_params.vacuum_length[0])
            && (x_cell[0]<prof_params.vacuum_length[0]+prof_params.length_params_x[0]) ) {
            fx = 1.0;
        } else {
            fx = 0.0;
        }
        
        // y-direction
        if (   (x_cell[1]>prof_params.vacuum_length[1])
            && (x_cell[1]<prof_params.vacuum_length[1]+prof_params.length_params_y[0]) ) {
            fy = 1.0;
        } else {
            fy = 0.0;
        }
        
        // x-y direction
        return fx*fy;
        
    }// constant
    
    
    // Harris velocity profile
    // -----------------------
    // Harris density profile: used for reconnection
    // ---------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[0] : characteristic length of the Harris profile
    // double_params[0]   : density nb parameter
    else if (prof_params.profile=="harris") {
        
        double nb = prof_params.double_params[0];
        double L  = prof_params.length_params_y[0];
        double y0 = prof_params.length_params_y[1];
        
        return (1.0 - pow(tanh((x_cell[1]-y0)/L),2)) / (nb + 1.0/pow(cosh((x_cell[1]-y0)/L),2));
    }

}

