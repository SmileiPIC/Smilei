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

//! this structure holds all the possible paraeters for phase diagnostics. Then every DiagnosticPhaseXXXXX will pick the ones that fit
struct phaseStructure {
	//!string defining the kind oh phase projections
    std::vector<std::string> kind;

    //! phase output every (every phase diagnostic must have this)
    unsigned int every;
	
	//! vector of pointer to species on which the phase diag will be applied (if omitted, it will be for all)
	std::vector<std::string> species;
	
    //! minimum position
	std::vector<double> pos_min;
    //! max position
	std::vector<double> pos_max;
    //! number of positions
	std::vector <unsigned int> pos_num;

    //! minimum momentum
	std::vector<double> mom_min;
    //! max momentum
	std::vector<double> mom_max;
    //! number of momenta
	std::vector <unsigned int> mom_num;
	
    //! minimum Lorentz factor
	std::vector<double> lor_min;
    //! max Lorentz factor
	std::vector<double> lor_max;
    //! number of Lorentz factors
	std::vector <unsigned int> lor_num;
	
};

//! this structure contains the definition of a probe0D
struct probe0DStructure {
    //! probe0D output every (every probe1D diagnostic must have this)
    unsigned int every;
    
    //! rearranged positions for the probes 0d
    std::vector<double> pos;

};

//! this structure contains the definition of a probe1D
struct probe1DStructure {
    //! probe1D output every (every probe1D diagnostic must have this)
    unsigned int every;
    //! start position of the 1D probe
    std::vector<double> posStart;
    //! end position of the 1D probe
    std::vector<double> posEnd;
    //! number of probes between posStart and posEnd
    unsigned int number;
};

//! this structure contains the definition of a probe1D
struct probe2DStructure {
    //! probe2D output every (every probe1D diagnostic must have this)
    unsigned int every;
    //! center position of the 2D probe
    std::vector<double> posCenter;
    //! end position of the first axec of the 2D probe
    std::vector<double> posEndFirst;
    //! end position of the second axec of the 2D probe
    std::vector<double> posEndSecond;
    //! number of probes between posCenter and posEndFirst and between posCenter and posEndSecond
    std::vector<unsigned int> number;
};

// ---------------------------------------------------------------------------------------------------------------------
//! DiagParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class DiagParams {

public:
    //! Creator for DiagParams
    DiagParams(InputData &, PicParams&);

    //! field dump output
    unsigned int fieldDump_every;
    
    //! time-averaged field dump output
    unsigned int avgfieldDump_every;
    
    //! number of time-steps for time-averaging of fields
    unsigned int ntime_step_avg;

    //! particle dump output
    unsigned int particleDump_every;

    //! scalar output every scalar_every (namelist group "diagnostic scalar" key "every")
    unsigned int scalar_every;

    //! vector of 0D probes
    std::vector<probe0DStructure> probe0DStruc;

    //! rearranged positions for the probes 0d
    std::vector<std::vector<double> > ps_0d_coord;

    //! vector of 1D probes
    std::vector<probe1DStructure> probe1DStruc;
    
    //! vector of 2D probes
    std::vector<probe2DStructure> probe2DStruc;
    
    //! every for the standard pic timeloop output
    unsigned int print_every;
	
	//! vector containing phase1D structures
	std::vector<phaseStructure> vecPhase;
};

#endif
