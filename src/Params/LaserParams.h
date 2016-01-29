/*! @file LaserParams.h

  @brief LaserParams.h is the class that hold the laser parameters and can read from a file the namelist

  @date 2013-02-15
*/

#ifndef LaserParams_H
#define LaserParams_H

#include "PyTools.h"

#include <vector>
#include <string>

class Params;

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Laser Profile
// ---------------------------------------------------------------------------------------------------------------------
struct LaserProfileStructure {
    
    //! Profile profile
    std::string profile;
    
    //! in case profile is give in Python
    PyObject *py_profile;
    
    //! int vector for profile parameters
    std::vector<int> int_params;
    
    //! double vector for profile parameters
    std::vector<double> double_params;
    
};


// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Laser
// ---------------------------------------------------------------------------------------------------------------------
struct LaserStructure {
    
    //! Side (west/east) from which the laser enters the box
    std::string boxSide;
    
    //! Laser field amplitude
    double a0;
    
    //! Laser angular frequency
    double omega0;
    
    //! Laser temporal chirp (for now assumed constant, maybe later use a fct?)
    double tchirp;
    
    //! Laser angle
    double angle;
    
    //! focus
    std::vector<double> focus;
    
    //! control parameter if laser is focused
    bool isFocused;
    
    //! Laser delta (ellipticity/polarization parameter)
    double delta;
    
    //! time profile 
    LaserProfileStructure profile_time;
    
    //! transverse profile
    LaserProfileStructure profile_transv;

    //! time-delay used when the laser as non-normal incidence
    double delay;
};



// ---------------------------------------------------------------------------------------------------------------------
//! LaserParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class LaserParams {

public:
    //! Creator for LaserParams
  LaserParams(Params&, bool);
    void print();
    void check_params( unsigned int n_laser );

    //! laser parameters
    std::vector<LaserStructure> laser_param;
    
    
};

#endif
