#include "Laser.h"

using namespace std;



Laser::Laser(double sim_time, std::vector<double> sim_length, LaserStructure laser_param) {

    pi_ov_2 = 0.5 * M_PI;

    laser_struct = laser_param;
    a0_delta_y_  = laser_struct.a0 * laser_struct.delta;
    a0_delta_z_  = laser_struct.a0 * sqrt(1.0-pow(laser_struct.delta,2));

    type_of_time_profile  = laser_struct.time_profile;
    type_of_y_profile  = laser_struct.y_profile;
    int_params            = laser_struct.int_params;
    double_params         = laser_struct.double_params;
    y_params         = laser_struct.y_params;

    // -------------------------------------------------------
    // Reconstruction of the vector int_params & double_params
    // and definition of the parameters
    // -------------------------------------------------------

    // CONSTANT time-profile
    // double_params[0]: duration of the pulse (default = sim_time)
    // double_params[1]: delay before the pulse is turned on (default = 0)
    if (type_of_time_profile=="constant") {

        if (double_params.size()<1) {
            double_params.resize(2);
            double_params[0] = sim_time;
        }
        else {
            double_params.resize(2);
        }
    }

    // SIN2 time-profile
    // double_params[0]: rise/fall-time of the sine-square profile (default = sim_time/2)
    // double_params[1]: duration of constant plateau at maximum (default = 0)
    // double_params[2]: delay before the pulse is turned on (default = 0)
    else if (type_of_time_profile=="sin2") {

        if (double_params.size()<1) {
            double_params.resize(3);
            double_params[0] = sim_time/2.0;
        }
        else {
            double_params.resize(3);
        }
    }
// gauss time-profile
// double_params[0]: Time at which maximum intensity is reached at the x=0 boundary.
// double_params[1]: Longitudinal FWHM of the pulse intensity.
    else if (type_of_time_profile=="gauss" ) {

        if (double_params.size()<1) {
            double_params.resize(2);
            double_params[0] = sim_time/2.0;
        }
        else {
            double_params.resize(2);
        }
    }

    else {
        ERROR("Laser profile " << type_of_time_profile <<  " not defined");
    }// ENDIF type_of_time_profile

//Constant transverse profile
   if (type_of_y_profile=="constant") {

        if (y_params.size()>0) {
            double_params.resize(0);
        }
    }
    else if (type_of_y_profile=="gauss"){

        if (double_params.size()<1) {
            double_params.resize(1);
            double_params[0] = sim_length[1]/4.0;
        }
        else {
            double_params.resize(1);
        }

    }
}


double Laser::time_profile(double time_dual) {


    // CONSTANT time-profile
    // double_params[0]: duration of the pulse (default = sim_time)
    // double_params[1]: delay before the pulse is turned on (default = 0)
    if (type_of_time_profile=="constant") {

        // delay before pulse
        if (time_dual<=double_params[1]) {
            return 0.0;
        }
        // plateau
        else if (time_dual<=double_params[1]+double_params[0]) {
            return 1.0;
        }
        // after the pulse
        else {
            return 0.0;
        }
    }


// SIN2 time-profile
// double_params[0]: rise/fall-time of the sine-square profile (default = sim_time/2)
// double_params[1]: duration of constant plateau at maximum (default = 0)
// double_params[2]: delay before the pulse is turned on (default = 0)
    else if (type_of_time_profile=="sin2") {

        // delay before pulse
        if (time_dual<=double_params[2]) {
            return 0.0;
        }
        // sin2 rise
        else if (time_dual<=double_params[2]+double_params[0]) {
            return pow( sin( pi_ov_2 * (time_dual-double_params[2]) / double_params[0] ) , 2 );
        }
        // plateau
        else if (time_dual<=double_params[2]+double_params[0]+double_params[1]) {
            return 1.0;
        }
        // sin2 fall
        else if (time_dual<=double_params[2]+double_params[0]+double_params[1]+double_params[0]) {
            return pow( cos(pi_ov_2 * (time_dual-(double_params[2]+double_params[0]+double_params[1]))
                            / double_params[0] ) , 2 );
        }
        // after the pulse
// gauss time-profile
// double_params[0]: Time at which maximum intensity is reached at the x=0 boundary.
// double_params[1]: Longitudinal FWHM of the pulse intensity.
    else if (type_of_time_profile=="gauss") {
        //exp(-2*log(2)*pow((time_dual-double_params[0]) / double_params[1] , 2)); //Gaussian longitudinal profil of the intensity
        //FWHM(Intensity) = FWHM(Field^2) = FWHM(Field)/sqrt(2) ==>
        return exp(log(2)*pow((time_dual-double_params[0]) / double_params[1] , 2)); //Gaussian longitudinal profil of the field as required
         
    } 
    else {
            return 0.0;
        }
    }


    else
        return 0.0;
}

double Laser::y_profile(double dfa) {

    if (type_of_y_profile=="constant") {
        return 1.0;
    }
    else if (type_of_y_profile=="gauss"){
        return exp(- pow(dfa / y_params[0] , 2));
    }
    else
        return 1.0;
}
