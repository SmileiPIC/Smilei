/*! @file DiagParams.h

  @brief DiagParams.h is the class that hold the simulation parameters and can read from a file the namelist

  @author tommaso vinci
  @date 2013-02-15
*/

#ifndef DIAGPARAMS_H
#define DIAGPARAMS_H

#include <vector>
#include <string>

#include "InputData.h"
#include "PicParams.h"

// ---------------------------------------------------------------------------------------------------------------------
//! DiagParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------

struct phaseStructure {
	//!string defining the kind oh phase
	std::string kind;

    //! phase 1D output every probe_every (namelist group "diagnostic phase1d" key "every")
    unsigned int every;
	
	//! vector of pointer to species on which the phase diag will be applied
	std::vector<std::string> species;
	
};

class DiagParams {

public:
    //! Creator for DiagParams
    DiagParams(InputData &, PicParams&);

    //! field dump output
    unsigned int fieldDump_every;

    //! particle dump output
    unsigned int particleDump_every;

    //! scalar output every scalar_every (namelist group "diagnostic scalar" key "every")
    unsigned int scalar_every;

    //! map output every map_every (namelist group "diagnostic map" key "every")
    unsigned int map_every;

    //! scalar output every probe_every (namelist group "diagnostic probe0d" key "every")
    unsigned int probe0d_every;
		
	//! number of 1D probes
    unsigned int n_probe1d;

    std::vector<std::vector<double> > ps_coord;
    std::vector<std::vector<std::vector<double> > > ps_1d_coord;
    std::vector<unsigned int> ps_1d_every;
    std::vector<unsigned int> ps_1d_res;

    //! every for the standard pic timeloop output
    unsigned int print_every;
	
	//! vector containing phase1D structures
	std::vector<phaseStructure> vecPhase;
};

#endif
