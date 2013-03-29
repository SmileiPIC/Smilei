#include "Laser.h"

Laser::Laser(LaserStructure laser_param) {
	laser_struct=laser_param;
	a0_delta_y_= laser_struct.a0 * laser_struct.delta;
	a0_delta_z_= laser_struct.a0 * sqrt(1.0-pow(laser_struct.delta,2));
}

double Laser::time_profile(double time_dual) {
	if (laser_struct.time_profile=="constant") {
		if (time_dual<=laser_struct.double_params[0]) {
			return 1.0;
		} else {
			return 0.0;
		}
        
	} else
        return 0.0;
}
