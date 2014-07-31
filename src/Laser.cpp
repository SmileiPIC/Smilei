#include "Laser.h"

using namespace std;



Laser::Laser(double sim_time, double sim_length_y, LaserStructure laser_param) {

    pi_ov_2 = 0.5 * M_PI;
    
    laser_struct = laser_param;
    a0_delta_y_  = laser_struct.a0 * laser_struct.delta;
    a0_delta_z_  = laser_struct.a0 * sqrt(1.0-pow(laser_struct.delta,2));

    type_of_time_profile   = laser_struct.time_profile;
    int_params             = laser_struct.int_params;
    double_params          = laser_struct.double_params;
    type_of_transv_profile = laser_struct.transv_profile;
    int_params_transv      = laser_struct.int_params_transv;
    double_params_transv   = laser_struct.double_params_transv;

    // -------------------------------------------------------
    // LASER TIME PROFILE INFOS
    // Reconstruction of the vector int_params & double_params
    // and definition of the laser time-profile parameters
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
    //double_params[1]: plateau
    //double_params[2]: delay
    //int_params[0]: gaussian cut-off
    else if(type_of_time_profile=="gaussian"){
        if(double_params.size()<1){
            double_params.resize(3);
            double_params[0]=sim_time/2.0;
            
        }
        else{
            double_params.resize(3);
        }
        if(int_params.size()<1){
            int_params.resize(1);
            int_params[0]=3;
        }
        else{
            int_params.resize(1);
        }
        
    }

    else {
        ERROR("Laser profile " << type_of_time_profile <<  " not defined");
    }// ENDIF type_of_time_profile

    
    // ---------------------------------------------------------------------
    // LASER TRANSVERSE PROFILE INFOS
    // Reconstruction of the vector int_params_transv & double_params_transv
    // and definition of the laser time-profile parameters
    // ---------------------------------------------------------------------
    
    // Plane-wave (default): no need of double_params_transv & int_params_transv
    if (type_of_transv_profile=="plane-wave") {
        MESSAGE(1,"Laser is a plane-wave");
    }
    
    // Gaussian or hyper-Gaussian transverse-profile
    // double_params_transv[0] : position in y of the center of the profile (default = length_sim[1]/2, middle of the box)
    // double_params_transv[1] : FWHM in intensity (default = length_sim[1]/4)
    // int_params_transv[0]    : order of the hyper-Gaussian profile  (default=2)
    else if (type_of_transv_profile=="gaussian") {
        MESSAGE(1,"Laser has a Gaussian or hyper-Gaussian transverse profile");
        if (double_params_transv.size()<2) {
            double_params_transv.resize(2);
            double_params_transv[0] = sim_length_y/2.0;
            double_params_transv[1] = sim_length_y/4.0;
        }
        if (int_params_transv.size()<1) {
            int_params_transv.resize(1);
            int_params_transv[0] = 2;
        }
    }
    
    // If transverse profile is not defined use the plane-wave as default
    else {
        type_of_transv_profile = "plane-wave";
        WARNING("Laser had no transverse profile defined: use plane-wave as default");
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
        else {
            return 0.0;
        }
    }
    //GAUSSIAN time-profile
    //double_params[0]: tau FWHM
    //double_params[1]: plateau
    //double_params[2]: delay
    //int_params[0]: gaussian cut-off
    else if(type_of_time_profile=="gaussian"){
        double fwhm=2*double_params[0];
        double sigma=fwhm/(2*sqrt(2*log(2)));
        double lt=int_params[0]*sigma;
        //delay before pulse
        if(time_dual<=double_params[2]){
            return 0.0;
        }
        //gaussian rise
        else if (time_dual<=double_params[2]+lt){
           
            return exp(-(time_dual-double_params[2]-lt)*(time_dual-double_params[2]-lt)/(2*sigma*sigma));

        }
        //plateau
        else if(time_dual<=double_params[2]+lt+double_params[1]){
            return 1.0;
        }
        //gaussian fall
        else if (time_dual<=double_params[2]+2*lt+double_params[1]){
            
            return exp(-(time_dual-double_params[2]-double_params[1]-lt)*(time_dual-double_params[2]-double_params[1]-lt)/(2*sigma*sigma));
            

        }

        //after pulse
        else{
            return 0.0;
        }
    }
    
    else
        return 0.0;
}//END laser::time_profile



double Laser::transverse_profile2D(double time_dual, double y) {
    
    
    // PLANE-WAVE
    if (type_of_transv_profile=="plane-wave") {
        return 1.0;
    }
    
    
    // GAUSSIAN or HYPER-GAUSSIAN PROFILE
    // double_params_transv[0] : position in y of the center of the profile
    // double_params_transv[1] : FWHM in intensity
    // int_params_transv[0]    : order of the hyper-Gaussian profile
    else if (type_of_transv_profile=="gaussian") {
        
        double sigma_N = pow(double_params_transv[1],int_params_transv[0]) / pow(2.0,int_params_transv[0]-1) / log(2.0);

        return exp( - pow(y-double_params_transv[0],int_params_transv[0]) / sigma_N);
    }
    
    else
        return 0.0;
}//END laser::transverse_profile2D