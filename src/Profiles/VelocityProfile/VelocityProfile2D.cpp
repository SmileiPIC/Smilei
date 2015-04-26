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
    // length_params_y[0] : position of first max in density
    // length_params_y[0] : position of second max in density
    // double_params[0]   : density nb parameter
    // double_params[1]   : maximum value of the magnetic field
    // double_params[2]   : electron over ion temperature ratio
    // int_params[0]      : 0 for ions, 1 for electrons
    else if (prof_params.profile=="harris") {
        
        if (prof_params.int_params.size()<1) {
            ERROR("For the Harris velocity profile 1 int_params has to be defined");
        }
        if (prof_params.double_params.size()<3) {
            ERROR("For the Harris velocity profile 3 double_params has to be defined");
        }
        if (prof_params.length_params_y.size()<3) {
            ERROR("For the Harris velocity profile 3 length_params_y have to be defined");
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
    

    // Harris velocity profile: used for reconnection
    // ---------------------------------------------
    // vacuum_length[0,1] : length of the vacuum region before the plasma in x & y directions (default is 0)
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[0] : characteristic length of the Harris profile
    // length_params_y[0] : characteristic length of the Harris profile
    // double_params[0]   : density nb parameter
    // double_params[1]   : maximum value of the magnetic field
    // double_params[2]   : electron over ion temperature ratio
    // int_params[0]      : 0 for ions, 1 for electrons
    
    else if (prof_params.profile=="harris") {
        
        double nb       = prof_params.double_params[0];
        double B0       = prof_params.double_params[1];
        double dB       = prof_params.double_params[2];
        double Te_ov_Ti = prof_params.double_params[3];
        int    N        = prof_params.int_params[0];
        double L   = prof_params.length_params_y[0];
        double y0  = prof_params.length_params_y[1];
        double y1  = prof_params.length_params_y[2];
        double sgl = prof_params.length_params_x[0];
        double x0  = prof_params.length_params_x[1];
        double x1  = prof_params.length_params_x[2];
        
        //std::cout << "nb" << nb << " sgl " << prof_params.length_params_x[0] << " " << sgl << std::endl;
        
        double Dx0 = (x_cell[0]-x0)/sgl;
        double Dx1 = (x_cell[0]-x1)/sgl;
        double Dy0 = (x_cell[1]-y0)/sgl;
        double Dy1 = (x_cell[1]-y1)/sgl;
        
        double Jz  = B0/L * pow(-Te_ov_Ti,N)/(1+Te_ov_Ti)
        *       ( pow(tanh((x_cell[1]-y0)/L),2) - pow(tanh((x_cell[1]-y1)/L),2) );
        
        //cout << dB << " " << sgl << " " << x0 << " " << x1 << endl;
        
        double dJz = 4.0*dB/sgl * (1.0-Dx0*Dx0-Dy0*Dy0) * exp(-Dx0*Dx0-Dy0*Dy0)
        -            4.0*dB/sgl * (1.0-Dx1*Dx1-Dy1*Dy1) * exp(-Dx1*Dx1-Dy1*Dy1);
        double n   = nb + pow(cosh((x_cell[1]-y0)/L),-2) + pow(cosh((x_cell[1]-y1)/L),-2);
        double v   = (Jz+dJz)/n;
        
        return v;
    }
    
    return 1;

}

