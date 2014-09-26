/*! @file LaserParams.h

  @brief LaserParams.h is the class that hold the laser parameters and can read from a file the namelist

  @author tommaso vinci
  @date 2013-02-15
*/

#ifndef LaserParams_H
#define LaserParams_H

#include <vector>
#include <string>

#include "InputData.h"
#include "PicParams.h"

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Laser
// ---------------------------------------------------------------------------------------------------------------------
struct LaserStructure {
    
    //! Side (west/east) from which the laser enters the box
    std::string boxSide;
    
    //! Laser field amplitude
    double a0;
    
    //! Laser angle
    double angle;
    
    //! focus
    std::vector<double> focus;
    
    //! Laser delta (ellipticity/polarization parameter)
    double delta;
    
    //! Laser profile
    std::string time_profile; //Longitudinal profile
    
    //! int vector for laser parameters
    std::vector<int> int_params;
    
    //! double vector for laser parameters
    std::vector<double> double_params; //Params for longitudinal profile
    
    //! Laser transverse profile
    std::string transv_profile;
    
    //! int vector for laser parameters
    std::vector<int> int_params_transv;
    
    //! double vector for laser parameters
    std::vector<double> double_params_transv;
};



// ---------------------------------------------------------------------------------------------------------------------
//! LaserParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class LaserParams {

public:
    //! Creator for LaserParams
    LaserParams(PicParams&, InputData &);

    //! initial number of laser pulses
    unsigned int n_laser;
    
    //! laser parameters
    std::vector<LaserStructure> laser_param;
    
    
};

#endif
