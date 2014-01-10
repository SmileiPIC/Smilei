
#ifndef LASER_H
#define LASER_H

#include "math.h"

#include "PicParams.h"

const double PI2 = 8.*atan(1.);

class Laser {
public:
	
	Laser(double, LaserStructure);
	
	LaserStructure laser_struct;

    double pi_ov_2;
    
	double a0_delta_y_;
	double a0_delta_z_;
	
    std::string         type_of_time_profile;
    std::vector<int>    int_params;
    std::vector<double> double_params;
    
	double time_profile(double);
	
	
};

#endif

