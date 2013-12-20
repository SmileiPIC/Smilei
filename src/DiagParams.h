/*! @file DiagParams.h 

  @brief DiagParams.h is the class that hold the simulation parameters and can read from a file the namelist

  @author tommaso vinci
  @date 2013-02-15
*/

#ifndef DIAGPARAMS_h
#define DIAGPARAMS_h

#include <vector>
#include <string>
#include "InputData.h"
#include "PicParams.h"

// ---------------------------------------------------------------------------------------------------------------------
//! DiagParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class DiagParams {
	
public:
	//! Creator for DiagParams
	DiagParams(InputData &, PicParams&);

	//! scalar output every scalar_every (namelist group "diagnostic scalar" key "every")
	unsigned int scalar_every;

	//! map output every map_every (namelist group "diagnostic map" key "every")
	unsigned int map_every;

	//! scalar output every probe_every (namelist group "diagnostic probe0d" key "every")
	unsigned int probe0d_every;
    std::vector<std::vector<double> > ps_coord;

};

#endif
