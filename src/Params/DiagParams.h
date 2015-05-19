/*! @file DiagParams.h

  @brief DiagParams.h is the class that hold the simulation parameters and can read from a file the namelist

  @date 2013-02-15
*/

#ifndef DIAGPARAMS_H
#define DIAGPARAMS_H

#include <vector>
#include <string>

#include "InputData.h"
#include "PicParams.h"
class Diagnostic;
class SmileiMPI;

//! this structure holds all the possible paraeters for phase diagnostics. Then every DiagnosticPhaseXXXXX will pick the ones that fit
struct phaseStructure {
	//!string defining the kind oh phase projections
    std::vector<std::string> kind;

    //! phase output every (every phase diagnostic must have this)
    unsigned int every;

    //! phase output from tmin
    double tmin;

    //! phase output to tmin
    double tmax;

    //! compression level using zlib [0-9] (0 deactvate compression)
    unsigned int deflate;

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


// ---------------------------------------------------------------------------------------------------------------------
//! DiagParams class: holds all the properties of the simulation that are read from the input file
// ---------------------------------------------------------------------------------------------------------------------
class DiagParams {

public:
    //! Creator for DiagParams
    DiagParams(Diagnostic&, PicParams&, InputData&, SmileiMPI *);

    void initScalars(Diagnostic&, PicParams&, InputData&);
    
    void initProbes(Diagnostic&, PicParams&, InputData&, SmileiMPI *);
    
    void initPhases(Diagnostic&, PicParams&, InputData&, SmileiMPI *);
    
    void initParticles(Diagnostic&, PicParams&, InputData&);
    
    //! field dump output
    unsigned int fieldDump_every;
    
    //! time-averaged field dump output
    unsigned int avgfieldDump_every;
    
    //! number of time-steps for time-averaging of fields
    unsigned int ntime_step_avg;

    //! particle dump output
    unsigned int particleDump_every;

    //! scalar output every scalar_every (namelist group "diag_scalar" key "every")
    unsigned int scalar_every;

    double scalar_tmin;
    double scalar_tmax;


    //! list of vars for scalars to be written (empty means all)
    std::vector<std::string> scalar_vars;

    //! scalar output precision
    unsigned int scalar_precision;
    
    //! every for the standard pic timeloop output
    unsigned int print_every;
	
	//! vector containing phase1D structures
	std::vector<phaseStructure> vecPhase;
	
	//! Method to find the numbers of requested species, sorted, and duplicates removed
    static std::vector<unsigned int> FindSpecies(std::vector<std::string>, PicParams&);

};

#endif
