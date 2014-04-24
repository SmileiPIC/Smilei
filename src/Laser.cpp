#include "Laser.h"

using namespace std;



Laser::Laser(double sim_time, LaserStructure laser_param) {

    pi_ov_2 = 0.5 * M_PI;

    laser_struct = laser_param;
    a0_delta_y_  = laser_struct.a0 * laser_struct.delta;
    a0_delta_z_  = laser_struct.a0 * sqrt(1.0-pow(laser_struct.delta,2));

    type_of_time_profile  = laser_struct.time_profile;
    int_params            = laser_struct.int_params;
    double_params         = laser_struct.double_params;

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
    //GAUSSIAN time-profile
    //double_params[0]: tau FWHM
    //double_params[1]: delay
    //int_params[0]: gaussian cut-off
    else if(type_of_time_profile=="gaussian"){
        if(double_params.size()<1){
            double_params.resize(2);
            double_params[0]=sim_time/2.0;
            double_params[1]=0.0;
        }
        else{
            double_params.resize(2);
        }
        if(int_params.size()<1){
            int_params.resize(1);
            int_params[0]=3.0;
        }
        else{
            int_params.resize(1);
        }
        
    }

    else {
        ERROR("Laser profile " << type_of_time_profile <<  " not defined");
    }// ENDIF type_of_time_profile


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
            DEBUG(pow( cos(pi_ov_2 * (time_dual-(double_params[2]+double_params[0]+double_params[1]))
                           / double_params[0] ) , 2 ))
            return pow( cos(pi_ov_2 * (time_dual-(double_params[2]+double_params[0]+double_params[1]))
                            / double_params[0] ) , 2 );
        }
        // after the pulse
        else {
            return 0.0;
        }
    }
    //GAUSSIAN time-profile
    //double_params[0]: tau FWHM
    //double_params[1]: delay
    //int_params[0]: gaussian cut-off
    else if(type_of_time_profile=="gaussian"){
        double fwhm=2*double_params[0];
        double sigma=fwhm/(2*sqrt(2*log(2)));
        double lt=int_params[0]*sigma;
        //delay before pulse
        if(time_dual<=double_params[1]){
            return 0.0;
        }
        //gaussian
        else if (time_dual<=double_params[1]+2*lt){
           
            double sol= exp(-(time_dual-double_params[1]-lt)*(time_dual-double_params[1]-lt)/(2*sigma*sigma));
//            DEBUG(sol)
            return sol;
        }
        //after pulse
        else{
            return 0.0;
        }
    }
    
    else
        return 0.0;
}
