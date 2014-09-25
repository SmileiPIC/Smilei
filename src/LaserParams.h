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
//! LaserParams class: holds all the properties of the lasers that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class LaserParams {

public:
    //! Creator for LaserParams
    LaserParams(PicParams&, InputData &);

};

#endif
