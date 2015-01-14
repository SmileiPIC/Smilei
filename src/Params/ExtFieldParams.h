/*! @file ExtFieldParams.h

  @brief ExtFieldParams.h is the class that hold the laser parameters and can read from a file the namelist

  @date 2015-01-13
*/

#ifndef ExtFieldParams_H
#define ExtFieldParams_H

#include <vector>
#include <string>

#include "InputData.h"
#include "PicParams.h"

class ExtFieldProfile;

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each ExtField
// ---------------------------------------------------------------------------------------------------------------------
struct ExtFieldStructure {
    
    //! ExtField profile
    std::string profile; 
    
    //! int vector for external field parameters
    std::vector<int> int_params;
    
    //! double vector for external field parameters
    std::vector<double> double_params;
    
    //! double vector for external field parameters (lengths: will be multiplied by 2pi)
    std::vector<double> length_params;
    
};



// ---------------------------------------------------------------------------------------------------------------------
//! ExtFieldParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class ExtFieldParams {

public:
    //! Creator for ExtFieldParams
    ExtFieldParams(PicParams&, InputData &);

    std::string geometry;
    
    //! external fields parameters string is the name of the field and
    std::map<std::string, std::vector<ExtFieldStructure> > structs;
    
};

#endif
