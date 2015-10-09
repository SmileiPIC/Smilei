#include "LaserProfile.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// LASER CONSTRUCTOR
// ---------------------------------------------------------------------------------------------------------------------
LaserProfile::LaserProfile( Params &params, LaserParams &laser_params, unsigned int n_laser) :
laser_struct(laser_params.laser_param[n_laser])
{
    
    MESSAGE(1,"Creating laser " << n_laser);
             
    double sim_time= params.sim_time;
    vector<double> sim_length= params.sim_length;
        
    omega0_                = laser_struct.omega0;
    tchirp_                = laser_struct.tchirp;
    
    a0_delta_y_            = laser_struct.a0 * laser_struct.delta;
    a0_delta_z_            = laser_struct.a0 * sqrt(1.0-pow(laser_struct.delta,2));

    type_of_time_profile   = laser_struct.profile_time.profile;
    int_params             = laser_struct.profile_time.int_params;
    double_params          = laser_struct.profile_time.double_params;

    type_of_transv_profile = laser_struct.profile_transv.profile;
    int_params_transv      = laser_struct.profile_transv.int_params;
    double_params_transv   = laser_struct.profile_transv.double_params;
    


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
    // double_params[0]: rise/fall-time of the sine-square profile (default = sim_time/2), also corresponds to FWHM
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
    // python function probably empty
    else if (type_of_time_profile=="python") {
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
        MESSAGE(2,"Laser is a plane-wave");
    }
    
    // Gaussian or hyper-Gaussian transverse-profile
    // double_params_transv[0] : position in y of the center of the profile (default = length_sim[1]/2, middle of the box)
    // double_params_transv[1] : FWHM in intensity (default = length_sim[1]/4)
    // int_params_transv[0]    : order of the hyper-Gaussian profile  (default=2)
    else if (type_of_transv_profile=="gaussian") {
        MESSAGE(2,"Laser has a Gaussian or hyper-Gaussian transverse profile");
        if (double_params_transv.size()<2) {
            double_params_transv.resize(2);
            double_params_transv[0] = sim_length[1]/2.0;
            double_params_transv[1] = sim_length[1]/4.0;
        }
        if (int_params_transv.size()<1) {
            int_params_transv.resize(1);
            int_params_transv[0] = 2;
        }
    }
    
    // focused: laser beam with either an arbitrary incident angle or focus
    else if (type_of_transv_profile=="focused") {
        MESSAGE(2,"Laser has a Gaussian transverse profile with arbitrary focus or incidence");
        if (double_params_transv.size()<1) {
            double_params_transv.resize(1);
            double_params_transv[0] = sim_length[1]/4.0;
            WARNING("Laser waist redefined to sim_length[1]/4 = " << double_params_transv[0]);
        }
    }
    // python function probably empty
    else if (type_of_transv_profile=="python") {
    }
    
    // If transverse profile is not defined use the plane-wave as default
    else {
        type_of_transv_profile = "plane-wave";
        WARNING("Laser had no transverse profile defined: use plane-wave as default");
    }

}



// ---------------------------------------------------------------------------------------------------------------------
// DEFINITION OF THE TIME_PROFILE
// ---------------------------------------------------------------------------------------------------------------------
double LaserProfile::time_profile(double time_dual) {


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
    // double_params[0]: rise/fall-time of the sine-square profile (default = sim_time/2), also corresponds to FWHM
    // double_params[1]: duration of constant plateau at maximum (default = 0)
    // double_params[2]: delay before the pulse is turned on (default = 0)
    else if (type_of_time_profile=="sin2") {
        
        double tp = double_params[1];
        double t0 = double_params[2];
        double T  = double_params[0];

        // delay before pulse
        if (time_dual<=t0) {
            return 0.0;
        }
        // sin2 rise
        else if (time_dual<=t0+T) {
            return abs( sin(M_PI_2*(time_dual-t0)/T) );
        }
        // plateau
        else if (time_dual<=t0+T+tp) {
            return 1.0;
        }
        // sin2 fall
        else if (time_dual<=t0+T+tp+T) {
            return abs( cos(M_PI_2*(time_dual-(t0+T+tp))/T) );
        }
        // after the pulse
        else {
            return 0.0;
        }
    }

    //GAUSSIAN time-profile
    //double_params[0]: tau FWHM : Field or intensity ? I guess intensity ?
    //double_params[1]: plateau : 
    //double_params[2]: delay : Time before the head of the laser pulse enters the box.
    //int_params[0]: gaussian cut-off : How far ahead from the center the pulse starts being non zero in number of sigma.
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
    
    // python function
    else if (type_of_time_profile=="python") {
        return PyTools::runPyFunction(laser_struct.profile_time.py_profile,time_dual);
    } 
    else {
        return 0.0;
    }

}//END laser::time_profile



// ---------------------------------------------------------------------------------------------------------------------
// DEFINITION OF THE TRANSVERSE PROFILE IN 2D
// ---------------------------------------------------------------------------------------------------------------------
double LaserProfile::transverse_profile2D(double time_dual, double y) {
    
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
    
    // GAUSSIAN PROFILE with ARBITRARY INCIDENCE ANGLE & FOCUSING
    // double_params_transv[0] : FWHM in intensity at the focal spot (beam waist)
    // a cut is defined here according to the integer parameter int_params_transv[0] (default equal to 4)
    else if (type_of_transv_profile=="focused") {
        double cut=int_params_transv[0];
        if ( (y>=-cut)&&(y<=cut) ) {
            return exp( -0.5*pow(y,2) );
        }
        else {
            return 0.0;
        }
    }
    // python function
    else if (type_of_transv_profile=="python") {
        return PyTools::runPyFunction(laser_struct.profile_transv.py_profile, time_dual,y);
    } 
    
    else
        return 0.0;
}//END LaserProfile::transverse_profile2D


void LaserProfile::disabled()
{
    a0_delta_y_ = 0.;
    a0_delta_z_ = 0.;
}


