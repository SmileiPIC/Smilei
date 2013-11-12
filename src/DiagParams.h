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

// ---------------------------------------------------------------------------------------------------------------------
//! DiagParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class DiagParams {
	
public:
	//! Creator for DiagParams
	DiagParams();

	//! \todo Comment these 3 stuffs
	void parseInputData(InputData &, PicParams&);

	int scalar_every;
	int map_every;
	
};

#endif
