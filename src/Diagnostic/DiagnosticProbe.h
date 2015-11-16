/*
-----------------------------------------------------------------------
PROBE DIAGNOSTICS                   - Mickael & Julien ?? - 2014
                                    - Modified by F Perez - 04/2015
-----------------------------------------------------------------------

The probe diagnostics are used to interpolate fields at other locations
than the PIC grid. These locations can be set as a 0-D, 1-D or 2-D grid,
not necessarily parallel to the PIC grid.

In the input (namelist) file, each diagnostic is provided as follows:

# PROBE DIAGNOSTICS - interpolate the fields on a N-D arbitrary grid
# ---------------------------------------------------------------------------------
# every        => an integer : number of time-steps between each output
# time_range   => two floats : min and max times to output (all times if omitted)
# number       => N floats : number of grid points in each dimension
# pos          => N floats : position of the reference point
# pos_first    => N floats : position of the first point
# pos_second   => N floats : position of the second point

where N is the number of dimensions of the probe.
The arguments `pos`,  `pos_first` and `pos_second` define the positions of the
"ends" or "corners" of the grid.
For creating a 1-D grid `pos_second` should be omitted.
For creating a 0-D grid `pos_first` should also be omitted.


>> Example: 0-D probe in 1-D simulation
diag_probe
    every = 1
    pos   = 1.2
end

>> Example: 1-D probe in 1-D simulation
diag_probe
    every = 1
    pos       = 1.2
    pos_first = 5.6
    number    = 100
end

>> Example: 1-D probe in 2-D simulation
diag_probe
    every = 1
    pos       = 1.2  4.
    pos_first = 5.6  4.
    number    = 100
end

>> Example: 2-D probe in 2-D simulation
diag_probe
    every = 1
    pos        = 0.    0.
    pos_first  = 10.   0.
    pos_second = 0.    10.
    number     = 100   100
end


*/


#ifndef DiagnosticProbe_H
#define DiagnosticProbe_H

#include <cmath>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <hdf5.h>

#include "Tools.h"

#include "Species.h"
#include "Interpolator.h"
#include "Particles.h"

class Params;
class SmileiMPI;
class ElectroMagn;
class Field2D;

//! this class holds the point probe
class DiagnosticProbe {
    
public:
    
    //! the creator
    DiagnosticProbe(Params& params, SmileiMPI *smpi);
    
    ~DiagnosticProbe();
    
    //! run all probes
    void run(unsigned int timestep, ElectroMagn* EMfields, Interpolator* interp);
    
    //! return name of the probe based on its number
    std::string probeName(int p);

    //! vector containing the timesteps at which calculate each probe
    std::vector<unsigned int> every;

    std::vector<double> tmin;
    std::vector<double> tmax;
    double dt;
    
    //! fake particles acting as probes
    std::vector<Particles> probeParticles;
    
    //! number of fake particles for each probe diagnostic
    std::vector<unsigned int> nPart_total;
    
    //! each probe will write in a buffer
    std::vector< Field2D* > probesArray;
    
    std::vector<int> probesStart;

    //! Number of fields to save
    std::vector<int> nFields;
    
    //! List of fields to save
    std::vector<std::vector<std::string>> fieldname;
    
    //! Indices in the output array where each field goes
    std::vector<std::vector<unsigned int>> fieldlocation;
    
protected:
    //! E local fields for the projector
    LocalFields Eloc_fields;
    //! B local fields for the projector
    LocalFields Bloc_fields;
    //! J local fields for the projector
    LocalFields Jloc_fields;
    //! Rho local field for the projector
    double Rloc_fields;
    
};
#endif
