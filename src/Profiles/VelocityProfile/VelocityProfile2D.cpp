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
    // ---------------------------------------------------------------------------------
    // Charles magnetic field profile for Liang simulations
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_y[0] = position of the maximum magnetic field
    // length_params_y[1] = Length of the magnetic gradient
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile == "magexpansion") {
        //if (prof_params.int_params.size()<2)
        //    ERROR("two int_params must be defined for Charles velocity profile" );
        if (prof_params.double_params.size()<4)
            ERROR("four double_params must be defined for Charles velocity profile" );
        if (prof_params.length_params_y.size()<2)
            ERROR("two length_params_y must be defined for Charles velocity profile" );
    }     
    
    // ---------------------------------------------------------------------------------
    // Blob magnetic field profile for Liang simulations
    // double_params[0]   = background density
    // double_params[1]   = relative density variation
    // double_params[2]   = Total electrons pressure
    // double_params[3]   = Background magnetic field
    // double_params[4]   = Maximum magnetic field
    // double_params[5]   = phase argument
    // length_params_x[0] = x position of the maximum of the B-field
    // length_params_x[1] = Length of the magnetic gradient
    // length_params_y[2] = y position of the maximum of the B-field
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile == "blob") {
        if (prof_params.double_params.size()<6)
            ERROR("four double_params must be defined for Charles velocity profile" );
        if (prof_params.length_params_x.size()<1)
            ERROR("two length_params_x must be defined for Charles velocity profile" );
        if (prof_params.length_params_y.size()<2)
            ERROR("two length_params_y must be defined for Charles velocity profile" );
    } 
    
    
    // ---------------------------------------------------------------------------------
    // Cosinus of x velocity profile
    // length_params_x[0] = wavelength
    // length_params_x[1] = phase of the argument
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile == "cos_x") {
        //if (prof_params.int_params.size()<2)
        //    ERROR("two int_params must be defined for Charles velocity profile" );
        if (prof_params.length_params_x.size()<2)
            ERROR("two length_params_x must be defined for Charles velocity profile" );
    } 
    else if (prof_params.profile=="python") {
        DEBUG("it's a python profile");
    }
    else {
        ERROR("Profile " << prof_params.profile << " not defined");
    }
    
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
    
    // ------------------------
    // Charles velocity profile
    // ------------------------
    // vacuum_length[0]  : not used here
    // double_params[0]   = background density
    // double_params[1]   = Total plasma pressure at infinity P0 = n0*(Te + Ti +...)
    // double_params[2]   = background magnetic field
    // double_params[3]   = Maximum magnetic field
    // length_params_y[0] = position of the maximum magnetic field
    // length_params_y[1] = Length of the magnetic gradient
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile=="magexpansion") {
        double n0    = prof_params.double_params[0];
        double P0    = prof_params.double_params[1];
        double B0    = prof_params.double_params[2];
        double Bmax  = prof_params.double_params[3];
        double x0    = prof_params.length_params_y[0];
        double L     = prof_params.length_params_y[1];
        double x     = x_cell[1]-x0;
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
        else {	double Bm = Bmax;
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
    
    
    // ---------------------------------------------------------------------------------
    // Blob magnetic field profile for Liang simulations
    // double_params[0]   = background density
    // double_params[1]   = relative density variation
    // double_params[2]   = Total electrons pressure
    // double_params[3]   = Background magnetic field
    // double_params[4]   = Maximum magnetic field
    // double_params[5]   = phase argument
    // length_params_x[0] = x position of the maximum of the B-field
    // length_params_x[1] = Length of the magnetic gradient
    // length_params_y[2] = y position of the maximum of the B-field
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile=="blob") {
        double n0    = prof_params.double_params[0];
        double dn    = prof_params.double_params[1];
        double P0    = prof_params.double_params[2];
        double B0    = prof_params.double_params[3];
        double Bmax  = prof_params.double_params[4];
        double theta0= prof_params.double_params[5];
        double v0= prof_params.double_params[6];
        double x0    = prof_params.length_params_x[0];
        double y0    = prof_params.length_params_y[0];
        double L     = prof_params.length_params_x[1];
        double  r    = sqrt(pow(x_cell[0]-x0,2) + pow(x_cell[1]-y0,2));
        double ne    = n0*( 1 - dn/(pow(cosh(r/L),2)));
        double theta = atan((x_cell[1]-y0)/(x_cell[0]-x0));
        double vd;
        if (r<5*L) { vd = v0;}
        else {vd=v0;}
        double v;
        if (Bmax == 0.) {
            double Bm = sqrt(pow(B0,2) + 2*P0)-B0;	
            double v     = -2*Bm/(L*ne)*tanh(r/L) /(pow(cosh(r/L),2))*cos(theta + theta0) +vd;
        	if (abs(v)>1.0) ERROR("Velocity profile exceeding c");
            return v;
        }
        else {	double Bm = Bmax;
            double v     = -2*Bm/(L*ne)*tanh(r/L) /(pow(cosh(r/L),2))*cos(theta + theta0)+vd;
        	if (abs(v)>1.0) ERROR("Velocity profile exceeding c");
            return v;
        }
    }
    
    
    // Cos profileover the x direction
    // ------------------------
    // vacuum_length[0]  : not used here
    // length_params_x[0] = wavelength
    // length_params_x[1] = phase of the argument
    // ---------------------------------------------------------------------------------
    else if (prof_params.profile=="cos_x") {
        double L    = prof_params.length_params_x[0];
        double x0 = prof_params.length_params_x[1];
        double x    = x_cell[0]-x0;
        double pi   = 4*atan(1.);
        return cos(2*pi/L*x);
	}
    
    else if (prof_params.profile=="python") {
        PyObject *pyresult = PyObject_CallFunction(prof_params.py_profile, const_cast<char *>("dd"), x_cell[0], x_cell[1]);
        return PyTools::get_py_result(pyresult);
    }
    
    return 1;
    
}

