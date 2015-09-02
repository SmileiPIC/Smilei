/*! @file ExtFieldParams.h

  @brief ExtFieldParams.h is the class that hold the laser parameters and can read from a file the namelist

  @date 2015-01-13
*/

#ifndef ExtFieldParams_H
#define ExtFieldParams_H

#include <vector>
#include <string>

#include "Params.h"

class ExtFieldProfile;

// ---------------------------------------------------------------------------------------------------------------------
//! This structure contains the properties of each ExtField
// ---------------------------------------------------------------------------------------------------------------------
struct ExtFieldStructure : ProfileStructure {
    //! fields to which apply the exeternal field
    std::vector<std::string> fields;
};

// ---------------------------------------------------------------------------------------------------------------------
//! ExtFieldParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class ExtFieldParams {

public:
    //! we copy this from Params
    std::string geometry;
    
    //! Creator for ExtFieldParams
    ExtFieldParams(Params&);

    //! external fields parameters the key string is the name of the field and the value is a vector of ExtFieldStructure
    std::vector<ExtFieldStructure> structs;
    
    //! Antenna
    std::vector<ExtFieldStructure> antennas;
    
};

#endif
