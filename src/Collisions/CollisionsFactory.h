/*! @file CollisionsFactory.h

  @brief CollisionsFactory.h contains functions for the creation and initialization of Collisions objects

  @date 2015-03
*/

#ifndef CollisionsParams_H
#define CollisionsParams_H

#include <vector>

#include "InputData.h"
#include "PicParams.h"
#include "Species.h"
#include "Collisions.h"

// ----------------------------------------------------------------------
//! CollisionsFactory class: Creates Collisions objects by reading through the input file
// ----------------------------------------------------------------------
class CollisionsFactory {

public:
    //! Method that creates a vector of Collisions objects: one for each group in the input file.
    static std::vector<Collisions*> create(PicParams&, InputData&, std::vector<Species*>&);

};

#endif
