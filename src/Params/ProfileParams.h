/*! @file ProfileParams.h

  @brief ProfileParams.h is the class that hold the laser parameters and can read from a file the namelist

  @date 2015-01-13
*/

#ifndef ProfileParams_H
#define ProfileParams_H

#include <Python.h>
#include <vector>
#include <string>

class PicParams;


// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each Profile
// ---------------------------------------------------------------------------------------------------------------------
struct ProfileStructure {

    //! Profile profile
    std::string profile; 
    
    //! in case profile is give in Python
    PyObject *py_profile;

    //! int vector for profile parameters
    std::vector<int> int_params;
    
    //! double vector for profile parameters
    std::vector<double> double_params;

    //! double vector for profile parameters (x lengths: will be multiplied by 2pi)
    std::vector<double> length_params_x;
    
    //! double vector for profile parameters (y lengths: will be multiplied by 2pi)
    std::vector<double> length_params_y;
    
    //! double vector for profile parameters (z lengths: will be multiplied by 2pi)
    std::vector<double> length_params_z;
    
};


// ---------------------------------------------------------------------------------------------------------------------
//! for the species we need an additional variable
// ---------------------------------------------------------------------------------------------------------------------
struct ProfileSpecies : ProfileStructure {
    //! vacuum lengths
    std::vector<double> vacuum_length;
};


// ---------------------------------------------------------------------------------------------------------------------
//! ProfileParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class ProfileParams {

public:
    //! Creator for ProfileParams
    ProfileParams(PicParams&);

    //! we copy this from picparams
    std::string geometry;
    
};




#endif
