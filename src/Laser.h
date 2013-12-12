
#ifndef LASER_H
#define LASER_H

#include "math.h"

#include "PicParams.h"

const double PI2 = 8.*atan(1.);

class Laser {
public:
	
	Laser(LaserStructure);
	
	LaserStructure laser_struct;

	double a0_delta_y_;
	double a0_delta_z_;
	
	double time_profile(double);
	
	
};

#endif

