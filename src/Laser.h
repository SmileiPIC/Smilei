#ifndef LASER_H
#define LASER_H

#include <cmath>

#include "PicParams.h"

const double PI2 = 8.*atan(1.);

class Laser {
public:

    Laser(double, std::vector<double>, LaserStructure);

    LaserStructure laser_struct;

    double pi_ov_2;

    double a0_delta_y_;
    double a0_delta_z_;

    std::string         type_of_time_profile;
    std::string         type_of_y_profile;
    std::vector<int>    int_params;
    std::vector<double> double_params;
    std::vector<double> y_params;

    double time_profile(double);
    double y_profile(double);


};

#endif

